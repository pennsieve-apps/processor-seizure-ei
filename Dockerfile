FROM --platform=linux/amd64 python:3.12.16-slim-bookworm

WORKDIR /processor

# Pinned Python dependencies (binary wheels only, no apt packages), plus the Lambda runtime
# interface client for running as a Lambda function.
COPY requirements.txt /processor/requirements.txt
RUN pip install --no-cache-dir --only-binary=:all: -r /processor/requirements.txt \
    && find /usr/local/lib/python3.12 -depth -type d \( -name __pycache__ -o -name tests \) -prune -exec rm -rf {} +

COPY processor/ /processor/processor/
COPY VERSION /processor/VERSION

ENV PYTHONPATH="/processor" \
    PYTHONUNBUFFERED=1 \
    PYTHONDONTWRITEBYTECODE=1 \
    MPLCONFIGDIR=/tmp/matplotlib \
    MNE_DONTWRITE_HOME=true \
    HOME=/tmp

# Runs as uid 1000 on the platform; nothing here needs root.
USER 1000

# Runs main() once as a container, or serves processor.handler as a Lambda function
# (see processor/entry.py).
ENTRYPOINT ["python3", "-m", "processor.entry"]
