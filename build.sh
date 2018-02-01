set -ex

# SET THE FOLLOWING VARIABLES
# docker hub username
USERNAME=mattjvincent
# image name
IMAGE=mouse_map_converter

docker build -t $USERNAME/$IMAGE:latest .



