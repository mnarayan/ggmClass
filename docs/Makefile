PYTHON = python
MARKDOWN = markdown2
PDFLATEX = pdflatex
BIBTEX = bibtex
MKDOCS = mkdocs
SVG2PDF = svg2pdf
NPM = npm
MATLAB = $(MATLABROOT)/bin/matlab
DOT = dot

site=doc/site/docs

mfiles=\
+standardize/standardize_cols.m\
+standardize/successive_norm.m

mfile_src=$(addprefix matlab/src/, $(mfiles))
mfile_doc=$(patsubst %.m,%.md,$(addprefix $(site)/mfiles/, $(mfiles)))
mds=$(wildcard doc/site/docs/*.md)

# Web components
node_modules=doc/site/.build/node_modules
node_modules-dir=$(node_modules)/.stamp
jquery=$(node_modules)/jquery/dist/jquery.min.js
fontawesome=$(node_modules)/font-awesome/css/font-awesome.min.css
bootstrap=$(node_modules)/bootstrap/dist/css/bootstrap.min.css


doc/figures/svg/%.pdf : doc/figures/svg/%.svg
	$(SVG2PDF) "$(<)" "$(@)"

$(site)/mfiles/%.md : matlab/%.m $(site)/mfiles/.stamp doc/matdoc.py doc/matdocparser.py
	mkdir -p $(dir $(@))
	$(PYTHON) doc/matdoc.py "$(<)" > "$(@)"

$(node_modules-dir):
	mkdir -p $(node_modules)
	touch -t 197301010000 $(node_modules-dir)

$(bootstrap): $(node_modules-dir)
	cd $(node_modules) ; $(NPM) install bootstrap@3.3.6

$(fontawesome): $(node_modules-dir)
	cd $(node_modules) ; $(NPM) install font-awesome@4.5.0

$(jquery): $(node_modules-dir)
	cd $(node_modules) ; $(NPM) install jquery@2.1.4

doc/site/site/index.html : doc/site/mkdocs.yml \
  $(mfile_doc) $(mds) \
  $(bootstrap) $(fontawesome) $(jquery)
	mkdir -p doc/site/theme/{css,js,fonts}
	cd doc/site/theme/css ; ln -sfv ../../.build/node_modules/bootstrap/dist/css/bootstrap.min.css{,.map} ./
	cd doc/site/theme/css ; ln -sfv ../../.build/node_modules/font-awesome/css/font-awesome.min.css ./
	cd doc/site/theme/fonts ; ln -sfv ../../.build/node_modules/font-awesome/fonts/fontawesome-webfont.{eot,svg,ttf,woff,woff2} ./
	cd doc/site/theme/js ; ln -sfv ../../.build/node_modules/bootstrap/dist/js/{bootstrap.min.js,npm.js} ./
	cd doc/site/theme/js ; ln -sfv ../../.build/node_modules/jquery/dist/jquery.min.{js,map} ./
	cd doc/site ; $(MKDOCS) build --clean

doc-clean:
	rm -f doc/matdocparser.pyc doc/matdoc.pyc
	rm -f $(svg_tgt)

doc-distclean:
	rm -f doc/matconvnet-manual.pdf
	rm -rf doc/site/site
	rm -rf doc/site/.build
	rm -f $(mfile_doc)

doc-info:
	@echo "mds=$(mds)"
	@echo "mfile_src=$(mfile_src)"
	@echo "mfile_doc=$(mfile_doc)"
	@echo "svg_tgt=$(svg_tgt)"

doc-serve: doc
	cd doc/site ; $(MKDOCS) serve

