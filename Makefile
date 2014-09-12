
dum:
	@echo "Nothing done: specify target!"
	@echo "find : find all mardkown files in this directory"
	@echo "html : convert all markdown files in this directory to html files"

.PHONY: find

# find markdown/html files
find:
	find . -regex ".*.\.\(md\)"

html:
	for f in `find . -name "*.md"`; do \
		o=`echo $$f | sed 's/\.md/\.html/'`;  \
		markdown2 $$f > $$o; \
	done;


