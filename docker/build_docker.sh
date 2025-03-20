#!/bin/bash

# Name of the Docker image
IMAGE_NAME="bioinfo-toolkit"
IMAGE_TAG="latest"

# Build the Docker image for amd64 architecture
docker build --platform linux/amd64 -t $IMAGE_NAME:$IMAGE_TAG .
