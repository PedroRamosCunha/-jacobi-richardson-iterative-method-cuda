TARGET=nomedoprograma
CXX=c++
LD=c++
OBJS=main.cpp
nomedoprograma:$(OBJS)
	$(LD) -o $(TARGET) $(OBJS)
install:
	@install nomedoprograma /usr/local/bin/nomedoprograma
clean:
	rm -rf *.o
