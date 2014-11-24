function timeSound = playSound(pahandle, s)

PsychPortAudio('FillBuffer', pahandle, s);
timeSound = PsychPortAudio('Start', pahandle);
