TARGET = run

VPATH = .

CXX = gcc #-m64

SOURCES := $(foreach dir,$(VPATH),$(wildcard $(dir)/*))
SRCS = $(filter %.c,$(SOURCES))

OBJS = $(SRCS:%.c=%.o)
DEPS = $(OBJS:.o=.d)

LD_FLAGS = -L/home/xiajie/mpich2-1/lib -lmpich -lpthread -lrt -lm

CPPFLAGS = -Wall -I/home/xiajie/mpich2-1/include

%.d: %.c
	set -e; rm -f $@; \
	$(CXX) -MM $< > $@.$$$$; \
	sed 's,\($*\)\.o[ :]*,\1.o $@ : ,g' < $@.$$$$ > $@; \
	rm -f $@.$$$$

include $(OBJS:.o=.d)

$(OBJS): %.o:  %.c
	$(CXX) -c $(CPPFLAGS) $< -o $@

$(TARGET): $(OBJS)
	$(CXX) -o $@ $(OBJS) $(CPPFLAGS) $(LD_FLAGS)

all: $(TARGET)

.PHONY: clean

clean:
	rm -f $(OBJS) $(TARGET) $(DEPS) $(TESTS)

