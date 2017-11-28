#ifndef _SYSINFOR_H_
#define _SYSINFOR_H_

#include "utilities.h"
#include <Windows.h>

struct Sysinfo{
	LARGE_INTEGER frequency;
	LARGE_INTEGER t1, t2;
	float fFrameTimeQP;
	float fFrameTime;
	float fStartTime, fps;
	float currentTime;
	int iTotalFrames;
	char cTitleinfo[MAX_PATH];

	Sysinfo()
		:fFrameTime(0.0f),fFrameTimeQP(0.0f),currentTime(0.0),
		fStartTime(0.0f),fps(0.0f), iTotalFrames(0)
	{}
};

Sysinfo *sysinfoglobals;

#endif