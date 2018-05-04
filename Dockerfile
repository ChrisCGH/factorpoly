# Use any image as your base image, or "scratch"
# Add fwatchdog binary via https://github.com/openfaas/faas/releases/
# Then set fprocess to the process you want to invoke per request - i.e. "cat" or "my_binary"

FROM gcc:7.2

RUN echo "Pulling watchdog binary from Github." \
    && curl -sSL https://github.com/openfaas/faas/releases/download/0.6.13/fwatchdog > /usr/bin/fwatchdog \
    && chmod +x /usr/bin/fwatchdog

# install dependent packages
RUN apt update
RUN apt install -y libgmp-dev google-perftools
RUN ln -s /usr/lib/libtcmalloc_minimal.so.4 /usr/lib/libtcmalloc_minimal.so

# build handler
WORKDIR /c++/src/handler
COPY . .
RUN make 

FROM gcc:7.2

RUN apt update
RUN apt install -y libtcmalloc-minimal4
RUN ln -s /usr/lib/libtcmalloc_minimal.so.4 /usr/lib/libtcmalloc_minimal.so

# Add non root user
RUN groupadd app && useradd -g app app
RUN mkdir -p /home/app
RUN chown app /home/app

WORKDIR /home/app

COPY --from=0 /c++/src/handler/factorpoly .
COPY --from=0 /usr/bin/fwatchdog .
 
USER app

ENV fprocess=./factorpoly

CMD ["./fwatchdog"]
