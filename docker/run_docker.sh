#!/bin/bash

# Name of the Docker image
IMAGE_NAME="bioinfo-toolkit"
IMAGE_TAG="latest"

# Run the Docker container
docker run -it --rm --platform linux/amd64 -v $(pwd):/data $IMAGE_NAME:$IMAGE_TAG
