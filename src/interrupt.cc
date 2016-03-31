#include "interrupt.h"

cInterrupt::cInterrupt() : interrupted(false), active(true) {
	pthread_create(&this->thread, 0, &cInterrupt::loop, this);
}

void* cInterrupt::loop(void *obj) {
	while (reinterpret_cast<cInterrupt *>(obj)->active) reinterpret_cast<cInterrupt *>(obj)->readEv();
}

cInterrupt::~cInterrupt() {
	this->active = false;
	pthread_join(thread, 0);
}

void cInterrupt::readEv() {
	if (this->kb.getKeyState(KEY_ESC)) {
		this->interrupted = true;
		this->active = false;
	}
}
