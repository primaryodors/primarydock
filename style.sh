# Just a little utility to prevent indentation mismatches etc.

astyle --style=allman --indent=spaces -L -xW -w -xw src/*.cpp
astyle --style=allman --indent=spaces -L -xW -w -xw src/*.h
astyle --style=allman --indent=spaces -L -xW -w -xw src/classes/*.cpp
astyle --style=allman --indent=spaces -L -xW -w -xw src/classes/*.h
