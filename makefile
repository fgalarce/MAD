# Detect operating system in Makefile. Author: He Tao. Date: 2015-05-30.
ifeq ($(OS),Windows_NT)
  OSFLAG +=WIN32
  ifeq ($(PROCESSOR_ARCHITECTURE),AMD64)
    OSFLAG +=AMD64
  endif
  ifeq ($(PROCESSOR_ARCHITECTURE),x86)
    OSFLAG +=IA32
  endif
else
  UNAME_S := $(shell uname -s)
  ifeq ($(UNAME_S),Linux)
    OSFLAG +=LINUX
  endif
  ifeq ($(UNAME_S),Darwin)
    OSFLAG +=OSX
  endif
  UNAME_P := $(shell uname -p)
endif

GREEN=\033[0;32m
RED=\033[0;31m
BLUE=\033[0;34m
SHELL_COLOR=\033[0m

all:
	@echo "$(GREEN)----------------- TOOLS -----------------$(SHELL_COLOR)"
	make --directory=tools
	@echo "$(GREEN)----------------- IO --------------------$(SHELL_COLOR)"
	make --directory=io
	@echo "$(GREEN)----------------- SIMPLEX ---------------$(SHELL_COLOR)"
	make --directory=simplex
	@echo "$(GREEN)----------------- FEM -------------------$(SHELL_COLOR)"
	make --directory=fem
	@echo "$(GREEN)----------------- LA --------------------$(SHELL_COLOR)"
	make --directory=la
	@echo "$(GREEN)----------------- CALCULUS --------------$(SHELL_COLOR)"
	make --directory=calculus
	@echo "$(GREEN)----------------- PDE -------------------$(SHELL_COLOR)"
	make --directory=pde
	@echo "$(GREEN)----------------- ASSIMILATION ----------$(SHELL_COLOR)"
	make --directory=assimilation
	@echo "$(GREEN)-----------------------------------------$(SHELL_COLOR)"

clean:
ifeq ($(UNAME_S), Linux)
	rm fem/*.so
	rm simplex/*.so
	rm la/*.so
	rm assimilation/*.so
	rm pde/*.so
	rm calculus/*.so
	rm io/*.so
	rm tools/*.so
endif

ifeq ($(UNAME_S), Darwin)
	rm fem/*.dylib
	rm simplex/*.dylib
	rm la/*.dylib
	rm assimilation/*.dylib
	rm pde/*.dylib
	rm calculus/*.dylib
	rm io/*.dylib
	rm tools/*.dylib
endif
