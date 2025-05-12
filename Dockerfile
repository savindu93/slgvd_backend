# Python runtime
FROM python:3.12-slim

# Environmenta variables
ENV PYTHONDONTWRITEBYTECODE=1
ENV PYTHONUNBUFFERED=1

# Work dir
WORKDIR /slgvd_backend

# Install dependencies
COPY requirements.txt .
RUN pip install --upgrade pip && pip install -r requirements.txt

# Copy project
COPY . .

# Expose Daphne's port
EXPOSE 8080

# Run ASGI server
CMD ["daphne", "-b", "0.0.0.0", "-p", "8080", "backend.asgi:application"]


