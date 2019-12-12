FROM r-base:3.6.1

RUN apt-get update && apt-get install -y \
	python3.7 \
	python3-pip
