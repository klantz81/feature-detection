#ifndef INTERRUPT_H
#define INTERRUPT_H

#include "keyboard.h"

class cInterrupt {
  private:
	cKeyboard kb;
	bool active;
	pthread_t thread;
	
  protected:
    
  public:
	bool interrupted;
	  
	cInterrupt();
	~cInterrupt();
	static void* loop(void* obj);
	void readEv();
};

#endif