# Python runtime
FROM python:3.12-slim

# Environmenta variables
ENV PYTHONDONTWRITEBYTECODE=1
ENV PYTHONUNBUFFERED=1

# Work dir
WORKDIR /slgvd_backend

# Install dependencies
COPY requirements.txt .
RUN apt-get update && \
apt-get install -y gcc default-libmysqlclient-dev pkg-config curl && \
rm -rf /var/lib/apt/lists/*

# Download and make the Cloud SQL Proxy executable (for DB connection)
RUN curl -o /usr/local/bin/cloud-sql-proxy https://storage.googleapis.com/cloud-sql-connectors/cloud-sql-proxy/v2.12.0/cloud-sql-proxy.linux.amd64 && \
chmod +x /usr/local/bin/cloud-sql-proxy

RUN pip install --upgrade pip && pip install -r requirements.txt

# Copy project
COPY . .

# Expose Daphne's port
EXPOSE 8080

# Run ASGI server
CMD ["sh", "-c", "/usr/local/bin/cloud-sql-proxy --private-ip --port 3306 ${DB_CONNECTION_NAME} & sleep 5 && daphne -b 0.0.0.0 -p 8080 backend.asgi:application"]

# CMD ["daphne", "-b", "0.0.0.0", "-p", "8080", "backend.asgi:application"]


