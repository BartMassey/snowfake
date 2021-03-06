# Copyright © 2012 Bart Massey
# [This program is licensed under the "MIT License"]
# Please see the file COPYING in the source
# distribution of this work for license terms.

C99C = gcc -std=gnu99
CC = gcc
CFLAGS = -g -O4 -Wall -finline-functions
LIBS = -lm

snowfake: snowfake.c
	$(C99C) $(CFLAGS) -o snowfake snowfake.c $(LIBS)

clean:
	-rm -rf snowfake
