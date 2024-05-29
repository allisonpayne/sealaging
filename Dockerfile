# Use the rocker/tidyverse image as the base image
FROM rocker/tidyverse:latest

# Update package lists and install any additional dependencies if needed
RUN apt-get update && apt-get install -y \
    wget \
    unzip \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Copy directory to image
COPY . .
# Copy mark executable to default bin
COPY cluster/mark /usr/local/bin/mark

# Make mark executable
RUN chmod +x /usr/local/bin/mark

# Install the RMark package
RUN R -e "install.packages('RMark', repos='https://cloud.r-project.org/')"

# Specify commands to run when the container starts, if needed
CMD Rscript cluster/sealmark.R
