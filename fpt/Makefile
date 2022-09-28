all:
	javac -d bin $$(find src -name  '*.java')

win:
	dir /s /B *.java > sources.txt
	javac -d bin @sources.txt

run:
	cd bin && java main.MCSP
