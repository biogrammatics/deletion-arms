// Deletion Arms Designer - Frontend JavaScript

document.addEventListener('DOMContentLoaded', function() {
    // Elements
    const form = document.getElementById('designForm');
    const dropZone = document.getElementById('dropZone');
    const fileInput = document.getElementById('fastaFile');
    const dropZoneContent = document.getElementById('dropZoneContent');
    const fileInfo = document.getElementById('fileInfo');
    const fileName = document.getElementById('fileName');
    const clearFileBtn = document.getElementById('clearFile');
    const sequencesText = document.getElementById('sequencesText');
    const loading = document.getElementById('loading');
    const results = document.getElementById('results');
    const error = document.getElementById('error');
    const errorMessage = document.getElementById('errorMessage');
    const submitBtn = document.getElementById('submitBtn');
    const resetBtn = document.getElementById('resetBtn');
    const optimizationSlider = document.getElementById('optimizationBalance');
    const costWeightValue = document.getElementById('costWeightValue');
    const lengthWeightValue = document.getElementById('lengthWeightValue');
    const downloadFastaBtn = document.getElementById('downloadFasta');
    const downloadReportBtn = document.getElementById('downloadReport');
    const toggleEnzymesBtn = document.getElementById('toggleEnzymes');
    const enzymeSection = document.getElementById('enzymeSection');
    const enzymeArrow = document.getElementById('enzymeArrow');
    const enzymeList = document.getElementById('enzymeList');
    const selectAllEnzymes = document.getElementById('selectAllEnzymes');
    const deselectAllEnzymes = document.getElementById('deselectAllEnzymes');

    let currentFile = null;
    let lastFormData = null;
    let allEnzymes = [];

    // Slider update - single slider controls both weights (sum to 1)
    optimizationSlider.addEventListener('input', () => {
        const lengthWeight = parseFloat(optimizationSlider.value);
        const costWeight = (1 - lengthWeight).toFixed(1);
        costWeightValue.textContent = costWeight;
        lengthWeightValue.textContent = lengthWeight.toFixed(1);
    });

    // Load enzymes from API
    async function loadEnzymes() {
        try {
            const response = await fetch('/api/enzymes');
            const data = await response.json();
            allEnzymes = data.enzymes;
            renderEnzymeList();
        } catch (err) {
            enzymeList.innerHTML = '<p class="text-sm text-red-500">Failed to load enzymes</p>';
        }
    }

    function renderEnzymeList() {
        enzymeList.innerHTML = allEnzymes.map(enzyme => `
            <label class="flex items-center gap-2 text-sm cursor-pointer hover:bg-gray-50 p-1 rounded">
                <input type="checkbox" name="enzyme" value="${enzyme.name}" class="enzyme-checkbox rounded" checked>
                <span class="font-mono">${enzyme.name}</span>
                <span class="text-xs text-gray-400">$${enzyme.cost_per_unit.toFixed(3)}</span>
            </label>
        `).join('');
    }

    // Toggle enzyme section
    toggleEnzymesBtn.addEventListener('click', () => {
        enzymeSection.classList.toggle('hidden');
        enzymeArrow.classList.toggle('rotate-90');
    });

    // Select/Deselect all enzymes
    selectAllEnzymes.addEventListener('click', () => {
        document.querySelectorAll('.enzyme-checkbox').forEach(cb => cb.checked = true);
    });
    deselectAllEnzymes.addEventListener('click', () => {
        document.querySelectorAll('.enzyme-checkbox').forEach(cb => cb.checked = false);
    });

    // Get selected enzymes as comma-separated string
    function getSelectedEnzymes() {
        const checked = document.querySelectorAll('.enzyme-checkbox:checked');
        if (checked.length === 0) {
            return null; // No enzymes selected - will cause error
        }
        // Always return the checked enzymes
        return Array.from(checked).map(cb => cb.value).join(',');
    }

    // Load enzymes on page load
    loadEnzymes();

    // File upload handling
    dropZone.addEventListener('click', () => fileInput.click());

    dropZone.addEventListener('dragover', (e) => {
        e.preventDefault();
        dropZone.classList.add('border-primary', 'bg-blue-50');
    });

    dropZone.addEventListener('dragleave', () => {
        dropZone.classList.remove('border-primary', 'bg-blue-50');
    });

    dropZone.addEventListener('drop', (e) => {
        e.preventDefault();
        dropZone.classList.remove('border-primary', 'bg-blue-50');
        const files = e.dataTransfer.files;
        if (files.length > 0) {
            handleFile(files[0]);
        }
    });

    fileInput.addEventListener('change', () => {
        if (fileInput.files.length > 0) {
            handleFile(fileInput.files[0]);
        }
    });

    clearFileBtn.addEventListener('click', (e) => {
        e.stopPropagation();
        clearFile();
    });

    function handleFile(file) {
        currentFile = file;
        fileName.textContent = file.name;
        dropZoneContent.classList.add('hidden');
        fileInfo.classList.remove('hidden');
        sequencesText.value = '';
        sequencesText.disabled = true;
    }

    function clearFile() {
        currentFile = null;
        fileInput.value = '';
        dropZoneContent.classList.remove('hidden');
        fileInfo.classList.add('hidden');
        sequencesText.disabled = false;
    }

    // Form submission
    form.addEventListener('submit', async (e) => {
        e.preventDefault();

        // Validate input
        if (!currentFile && !sequencesText.value.trim()) {
            showError('Please upload a FASTA file or paste sequences.');
            return;
        }

        // Build form data
        const formData = new FormData();
        if (currentFile) {
            formData.append('file', currentFile);
        } else {
            formData.append('sequences', sequencesText.value);
        }
        formData.append('arm_length', document.getElementById('armLength').value);
        formData.append('half_site_min', document.getElementById('halfSiteMin').value);
        formData.append('half_site_max', document.getElementById('halfSiteMax').value);
        formData.append('max_designs', document.getElementById('maxDesigns').value);
        const lengthWeight = parseFloat(optimizationSlider.value);
        const costWeight = 1 - lengthWeight;
        formData.append('cost_weight', costWeight.toFixed(1));
        formData.append('length_weight', lengthWeight.toFixed(1));

        // Add selected enzymes (always send the checked list)
        const selectedEnzymes = getSelectedEnzymes();
        if (selectedEnzymes) {
            formData.append('available_enzymes', selectedEnzymes);
        } else {
            showError('Please select at least one enzyme.');
            return;
        }

        lastFormData = formData;

        // Show loading
        hideError();
        results.classList.add('hidden');
        loading.classList.remove('hidden');
        submitBtn.disabled = true;

        try {
            const response = await fetch('/api/design', {
                method: 'POST',
                body: formData
            });

            const data = await response.json();

            if (!response.ok) {
                throw new Error(data.detail || 'An error occurred');
            }

            displayResults(data);
        } catch (err) {
            showError(err.message);
        } finally {
            loading.classList.add('hidden');
            submitBtn.disabled = false;
        }
    });

    // Reset form
    resetBtn.addEventListener('click', () => {
        form.reset();
        clearFile();
        results.classList.add('hidden');
        hideError();
        optimizationSlider.value = '0.5';
        costWeightValue.textContent = '0.5';
        lengthWeightValue.textContent = '0.5';
        // Check all enzymes (default state)
        document.querySelectorAll('.enzyme-checkbox').forEach(cb => cb.checked = true);
    });

    // Download handlers
    downloadFastaBtn.addEventListener('click', () => downloadFile('/api/design/fasta', 'constructs.fa'));
    downloadReportBtn.addEventListener('click', () => downloadFile('/api/design/report', 'constructs_report.txt'));

    async function downloadFile(endpoint, filename) {
        if (!lastFormData) return;

        try {
            const response = await fetch(endpoint, {
                method: 'POST',
                body: lastFormData
            });

            if (!response.ok) throw new Error('Download failed');

            const blob = await response.blob();
            const url = window.URL.createObjectURL(blob);
            const a = document.createElement('a');
            a.href = url;
            a.download = filename;
            document.body.appendChild(a);
            a.click();
            window.URL.revokeObjectURL(url);
            a.remove();
        } catch (err) {
            showError('Download failed: ' + err.message);
        }
    }

    function showError(message) {
        errorMessage.textContent = message;
        error.classList.remove('hidden');
    }

    function hideError() {
        error.classList.add('hidden');
    }

    function displayResults(data) {
        results.classList.remove('hidden');

        // Check for no designs
        if (data.summary.total_designs === 0) {
            const summary = document.getElementById('summary');
            summary.innerHTML = `
                <div class="col-span-full bg-yellow-50 border border-yellow-200 rounded-lg p-6">
                    <h3 class="text-lg font-bold text-yellow-800 mb-2">No Valid Designs Found</h3>
                    <p class="text-yellow-700 mb-4">No construct designs could be generated with the current parameters. Try the following:</p>
                    <ul class="list-disc list-inside text-yellow-700 space-y-1">
                        <li><strong>Select more enzymes</strong> - Expand the enzyme selection to include more options</li>
                        <li><strong>Widen the half-site range</strong> - Increase the max or decrease the min distance</li>
                        <li><strong>Adjust arm length</strong> - Different arm lengths may have compatible half-sites</li>
                        <li><strong>Check your sequence</strong> - Ensure sufficient flanking regions around the deletion</li>
                    </ul>
                </div>
            `;
            document.getElementById('geneResults').innerHTML = '';
            return;
        }

        // Summary
        const summary = document.getElementById('summary');
        summary.innerHTML = `
            <div class="bg-blue-50 p-4 rounded-lg">
                <p class="text-2xl font-bold text-blue-700">${data.summary.total_genes}</p>
                <p class="text-sm text-blue-600">Genes</p>
            </div>
            <div class="bg-green-50 p-4 rounded-lg">
                <p class="text-2xl font-bold text-green-700">${data.summary.total_designs}</p>
                <p class="text-sm text-green-600">Total Designs</p>
            </div>
            <div class="bg-purple-50 p-4 rounded-lg">
                <p class="text-2xl font-bold text-purple-700">${data.summary.parameters.arm_length}</p>
                <p class="text-sm text-purple-600">Arm Length (bp)</p>
            </div>
            <div class="bg-orange-50 p-4 rounded-lg">
                <p class="text-2xl font-bold text-orange-700">${data.summary.parameters.cost_weight}/${data.summary.parameters.length_weight}</p>
                <p class="text-sm text-orange-600">Cost/Length Weight</p>
            </div>
        `;

        // Gene results
        const geneResults = document.getElementById('geneResults');
        geneResults.innerHTML = data.genes.map(gene => `
            <div class="bg-white rounded-lg shadow-md mb-6 overflow-hidden">
                <div class="bg-gray-800 text-white px-6 py-4">
                    <h3 class="text-lg font-bold">${gene.gene_name}</h3>
                    <p class="text-sm text-gray-300">Deletion size: ${gene.deletion_size} bp | ${gene.designs.length} design(s)</p>
                </div>
                <div class="overflow-x-auto">
                    <table class="w-full">
                        <thead class="bg-gray-50">
                            <tr>
                                <th class="px-4 py-3 text-left text-xs font-medium text-gray-500 uppercase">#</th>
                                <th class="px-4 py-3 text-left text-xs font-medium text-gray-500 uppercase">Enzyme(s)</th>
                                <th class="px-4 py-3 text-left text-xs font-medium text-gray-500 uppercase">Type</th>
                                <th class="px-4 py-3 text-left text-xs font-medium text-gray-500 uppercase">Down Arm</th>
                                <th class="px-4 py-3 text-left text-xs font-medium text-gray-500 uppercase">Up Arm</th>
                                <th class="px-4 py-3 text-left text-xs font-medium text-gray-500 uppercase">Total</th>
                                <th class="px-4 py-3 text-left text-xs font-medium text-gray-500 uppercase">Cost/unit</th>
                                <th class="px-4 py-3 text-left text-xs font-medium text-gray-500 uppercase">Details</th>
                            </tr>
                        </thead>
                        <tbody class="divide-y divide-gray-200">
                            ${gene.designs.map(design => `
                                <tr class="hover:bg-gray-50">
                                    <td class="px-4 py-3 text-sm">${design.design_num}</td>
                                    <td class="px-4 py-3 text-sm font-mono">
                                        ${design.single_enzyme
                                            ? design.enzyme1_name
                                            : `${design.enzyme1_name} + ${design.enzyme2_name}`}
                                    </td>
                                    <td class="px-4 py-3 text-sm">
                                        ${design.single_enzyme
                                            ? '<span class="bg-green-100 text-green-800 px-2 py-1 rounded text-xs">Single</span>'
                                            : '<span class="bg-yellow-100 text-yellow-800 px-2 py-1 rounded text-xs">Two</span>'}
                                    </td>
                                    <td class="px-4 py-3 text-sm">${design.downstream_arm_length} bp</td>
                                    <td class="px-4 py-3 text-sm">${design.upstream_arm_length} bp</td>
                                    <td class="px-4 py-3 text-sm font-medium">${design.total_arm_length} bp</td>
                                    <td class="px-4 py-3 text-sm">$${design.total_enzyme_cost.toFixed(4)}</td>
                                    <td class="px-4 py-3 text-sm">
                                        <button onclick="toggleDetails('${gene.gene_name}-${design.design_num}')"
                                                class="text-primary hover:underline text-xs">
                                            Show
                                        </button>
                                    </td>
                                </tr>
                                <tr id="details-${gene.gene_name}-${design.design_num}" class="hidden bg-gray-50">
                                    <td colspan="8" class="px-4 py-4">
                                        <div class="grid grid-cols-2 gap-4 text-sm">
                                            <div>
                                                <p class="font-medium text-gray-700">Enzyme 1: ${design.enzyme1_name}</p>
                                                <p class="text-gray-500">Recognition: ${design.enzyme1_recognition}</p>
                                                <p class="text-gray-500">Cost: $${design.enzyme1_cost_per_unit.toFixed(4)}/unit</p>
                                                <p class="text-gray-500">Position: ${design.downstream_half_site_position}</p>
                                            </div>
                                            <div>
                                                <p class="font-medium text-gray-700">Enzyme 2: ${design.enzyme2_name}</p>
                                                <p class="text-gray-500">Recognition: ${design.enzyme2_recognition}</p>
                                                <p class="text-gray-500">Cost: $${design.enzyme2_cost_per_unit.toFixed(4)}/unit</p>
                                                <p class="text-gray-500">Position: ${design.upstream_half_site_position}</p>
                                            </div>
                                        </div>
                                        ${design.stuffer_sequence ? `
                                            <div class="mt-4">
                                                <p class="font-medium text-gray-700">Stuffer Sequence (${design.stuffer_sequence.length} bp):</p>
                                                <p class="font-mono text-xs bg-gray-100 p-2 rounded mt-1 break-all">${design.stuffer_sequence}</p>
                                            </div>
                                        ` : ''}
                                        <div class="mt-4">
                                            <p class="font-medium text-gray-700">Construct Sequence (${design.construct_sequence.length} bp):</p>
                                            <p class="font-mono text-xs bg-gray-100 p-2 rounded mt-1 break-all max-h-32 overflow-y-auto">${design.construct_sequence}</p>
                                        </div>
                                    </td>
                                </tr>
                            `).join('')}
                        </tbody>
                    </table>
                </div>
            </div>
        `).join('');
    }

    // Global function for toggling details
    window.toggleDetails = function(id) {
        const details = document.getElementById('details-' + id);
        details.classList.toggle('hidden');
    };
});
