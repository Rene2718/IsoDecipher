FROM python:3.10-slim

# 1. Install system dependency
RUN apt-get update && apt-get install -y \
    curl \
    gnupg \
    git \
    procps \
    && curl https://sdk.cloud.google.com | bash

# 2. Define path
ENV PATH=$PATH:/root/google-cloud-sdk/bin

WORKDIR /app

# 3. Package installization
COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

# 4. Copy
COPY . .

# 5. Defaut the command
CMD ["python", "IsoDecipher/scripts/build_panel_features.py", "--help"]