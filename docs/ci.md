# CI Docker Image

The CI workflow uses a prebuilt Docker image hosted on GitHub Container Registry.

## Build and Push

```bash
# build the image
docker build -t ghcr.io/liusername/beidou-ci:latest -f docker/Dockerfile .

# login to GHCR (requires a token with write:packages)
echo "$GHCR_TOKEN" | docker login ghcr.io -u liusername --password-stdin

# push the image
docker push ghcr.io/liusername/beidou-ci:latest
```
