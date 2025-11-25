# Deletion Arms Designer - Docker Image
FROM python:3.11-slim

# Set working directory
WORKDIR /app

# Set environment variables
ENV PYTHONDONTWRITEBYTECODE=1
ENV PYTHONUNBUFFERED=1

# Install dependencies
COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

# Copy application code
COPY deletion_arms_designer.py .
COPY enzymes.tsv .
COPY web/ ./web/

# Expose port
EXPOSE 8000

# Run the application
CMD ["uvicorn", "web.app:app", "--host", "0.0.0.0", "--port", "8000"]
