SRCS := ${wildcard $(CURDIR)/clusterMon/*.cpp}
OBJS := $(SRCS:.cpp=.o)

PROG := $(CURDIR)/clusterMon/clusterMon
CXX := g++
CFLAGS := -Wall -O2 -std=c++1z

$(PROG): $(OBJS)
	$(CXX) -o $@ $^

%.o: %.cpp
	$(CXX) $(CFLAGS) -o $@ -c $<

clean:
	rm -f $(OBJS)
	rm -f $(PROG)
