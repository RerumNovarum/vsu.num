#!/bin/bash
cmake -DCMAKE_INSTALL_PREFIX=~/.local ..
make
make test
make install
