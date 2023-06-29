SRCS := ${wildcard *.cpp}
OBJS := $(SRCS:.cpp=.o)

PROG := hierarchical_lsh.exe
CXX := g++
CFLAGS := -Wall -O2 -std=c++1z -g

$(PROG): $(OBJS)
	$(CXX) -o $@ $^

%.o: %.cpp
	$(CXX) $(CFLAGS) -o $@ -c $<

clean:
	rm -f $(OBJS)
	rm -f $(PROG)
