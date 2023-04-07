SRCS := ${wildcard ./groundTruth/*.cpp}
OBJS := $(SRCS:.cpp=.o)

PROG := ./groundTruth/groundTruth
CXX := g++
CFLAGS := -Wall -O2 -std=c++1z

$(PROG): $(OBJS)
	$(CXX) -o $@ $^

%.o: %.cpp
	$(CXX) $(CFLAGS) -o $@ -c $<

clean:
	rm -f $(OBJS)
	rm -f $(PROG)
