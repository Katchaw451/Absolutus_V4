from absolutus_v4.fetch import human_gb, https_url


def test_https_url():
    assert https_url("s3://cellxgene-census-public-us-west-2/cell-census/2025-11-08/h5ads/x.h5ad", "us-west-2") == \
        "https://cellxgene-census-public-us-west-2.s3.us-west-2.amazonaws.com/cell-census/2025-11-08/h5ads/x.h5ad"
    assert https_url("https://example.org/a.h5ad") == "https://example.org/a.h5ad"


def test_human_gb():
    assert human_gb(None) == "unknown" and human_gb(2_500_000_000) == "2.50 GB"
