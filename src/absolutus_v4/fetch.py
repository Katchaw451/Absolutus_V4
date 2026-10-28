"""Small pure helpers for locating and sizing source .h5ad files."""
from __future__ import annotations


def https_url(uri: str, region: str = "us-west-2") -> str:
    """s3://bucket/key -> public virtual-hosted HTTPS URL (the Census bucket is public)."""
    if not uri.startswith("s3://"):
        return uri
    bucket, _, key = uri[5:].partition("/")
    return f"https://{bucket}.s3.{region}.amazonaws.com/{key}"


def human_gb(n_bytes) -> str:
    return "unknown" if n_bytes is None else f"{n_bytes / 1e9:.2f} GB"
