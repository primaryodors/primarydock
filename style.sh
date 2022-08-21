# Just a little utility to prevent indentation mismatches etc.

astyle -n --style=allman --indent=spaces -L -xW -w -xw src/*.cpp
# astyle -n --style=allman --indent=spaces -L -xW -w -xw src/*.h
astyle -n --style=allman --indent=spaces -L -xW -w -xw src/classes/*.cpp
astyle -n --style=allman --indent=spaces -L -xW -w -xw src/classes/*.h
