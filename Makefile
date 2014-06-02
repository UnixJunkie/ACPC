.PHONY: tags clean all conf build

all:
	obuild build

conf:
	obuild configure

build:
	obuild build

tags:
	otags `find . -name '*.ml' | egrep -v "myocamlbuild.ml|setup.ml"` -o TAGS

clean:
	obuild clean
	\rm -rf dist _build
