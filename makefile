run: script.mdl lex.py main.py matrix.py mdl.py display.py draw.py gmath.py yacc.py
	python main.py script.mdl

clean:
	rm *pyc *out parsetab.py

clear:
	rm -rf *pyc *out parsetab.py *ppm *png *gif anim/*png anim/*ppm
