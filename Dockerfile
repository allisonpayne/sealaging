# Use the rocker/r-ver image as the base image
FROM rocker/r-ver:latest

# Update package lists and install any additional dependencies if needed
RUN apt-get update && apt-get install -y \
    wget \
    unzip \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Set the working directory inside the container
WORKDIR /usr/src/app

# Download the mark.64.zip file
RUN wget http://www.phidot.org/software/mark/downloads/files/mark.64.zip

# Unzip the downloaded file
RUN unzip mark.64.zip

# Copy the contents to the default bin location and make them executable
RUN cp mark.64/mark64 /usr/local/bin/mark && \
    chmod +x /usr/local/bin/mark

# Install the RMark package
RUN R -e "install.packages('RMark', repos='https://cloud.r-project.org/')"

# Specify commands to run when the container starts, if needed
CMD [mark]