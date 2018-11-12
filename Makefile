init:
	pip install $(REQUIREMENTS)
	sed -i "1s/.*/$(path)python/" dgcg.py
	sed -i "4s/.*/python \"\$$BINPATH\/..\/dgcg\/dgcg.py\" $$\@/" ../bin/dgcg

#	pip install -r requirements.txt



REQUIREMENTS := astropy numpy scipy timeit 

path= \#!\/usr\/bin\/

python:
	sed -i "1s/.*/$(path)python/" dgcg.py 
	sed -i "4s/.*/python \"\$$BINPATH\/..\/dgcg\/dgcg.py\" $$\@/" ../bin/dgcg


python35:
	sed -i "1s/.*/$(path)python3.5/" dgcg.py 
	sed -i "4s/.*/python3.5 \"\$$BINPATH\/..\/dgcg\/dgcg.py\" $$\@/" ../bin/dgcg

python3:
	sed -i "1s/.*/$(path)python3/" dgcg.py
	sed -i "4s/.*/python3 \"\$$BINPATH\/..\/dgcg\/dgcg.py\" $$\@/" ../bin/dgcg



