# notebooks

[`contribute.ipynb`](contribute.ipynb) — one-click Colab runner for the dwell-time
screen. Sets a GPU runtime, installs the pip-only MD stack, and runs
[`screen/contribute_client.py`](../screen/contribute_client.py) against the public
contributor API. Linked from the screen's recruitment page (the "Open in Colab"
badge) so a volunteer with a Google account can donate GPU time without a local
install. The notebook is the *client*; the API it speaks is the contract — see
[`screen/`](../screen/) and the OpenAPI reference the site serves at `/screen/api`.
