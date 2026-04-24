#!/bin/bash

# Check if the user provided a filename
if [ $# -ne 1 ]; then
    echo "Usage: $0 <file_with_image_names>"
    exit 1
fi

IMAGE_LIST="$1"

# Check if the file exists
if [ ! -f "$IMAGE_LIST" ]; then
    echo "Error: File '$IMAGE_LIST' not found!"
    exit 1
fi

# Read each line and pull the corresponding Docker image
while IFS= read -r IMAGE; do
    if [ -n "$IMAGE" ]; then
        echo "Pulling image: $IMAGE"
        sudo docker pull "$IMAGE"
    fi
done < "$IMAGE_LIST"

echo "All images pulled successfully."
