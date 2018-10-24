(
// vars for NRT output
var outputPath, headerFormat, sampleFormat, numOutputChannels, sampleRate;

// bussing
var group, channels, garbage;

// imaging functions
var sinCosPanLaw, rotateMatrix, spatFilt, balanceMatrix, widthMatrix;

// filter functions
var sinc, lopass, hipass, bapass, bastop, loshelf, hishelf, dampFacFunc, sosFiltCoeffsFunc, cascadeCoeffsFunc;
var schroederCoeffFunc, regaliaMitraCoeff, combDelayFunc;

// synth def helper functions
var sampConvert, allPoleCoeffFunc, combT60Func, sampsToTime, equalPowMix, vibSinDelTimeFunc;

// synth defs
var playBufSynthDef, dReverb, feedbackEcho, slapbackEcho, sgsSynth, agsSynth,
sosHPFilter, sosLPFilter, sosBPFilter, sosBSFilter, firFilter, crossFade, convBasic, convEnv, universalComb, rotateImage, balanceImage, widthImage, vibTrem, combSpreader, rotarySpreader, weaverSSB, absShaper, sinShaper, cosShaper, disShaper;

// buffers
var buffersD, buffersC, sndPath, allSounds;

// note-calling
var score;
var buffer, bufL, bufR, start, dur;
var allPassBuffer, basicAmpEnv;
var loops;
var hiPass;

// note-calling functions
var getLeft, getRight;

/////////////// DEFINE NRT OUTPUT VARIABLES ///////////////

// set the NRT vars here...
outputPath = "Gaimari_final.wav".resolveRelative; // output file path
headerFormat = "WAV";                 // soundfile header format
sampleFormat = "int24";               // soundfile sample format
sampleRate = 44100;                   // sample rate
numOutputChannels = 2;                // stereo --> 2 channels

// create a score
score = CtkScore.new;

/////////////// DEFINE BUSSING VARIABLES ///////////////

// create busses in 2D array of [channel][stage]
// currently 5 channels each with 10 processing stages
channels = Array.fill2D(5, 10, {arg i; CtkAudio.new(2)});

garbage = CtkAudio.new(2);
// ... which is what we'll be sending!
// create a node group
group = CtkGroup.new;

///////////////    DEFINE IMAGING FUNCTIONS    ///////////////

// a function to return coefficients for Rotate- Matrix form
// position argument in degrees
rotateMatrix = { arg angleInDegrees = 0;
	var angleInRadians;
	var theta;

	angleInRadians = angleInDegrees/180*pi;

	theta = angleInRadians;

	Array2D.fromArray(
		2, // rows, outputs
		2, // columns, inputs
		[
			theta.cos, theta.sin,
			-1 * theta.sin, theta.cos
		]
	)
};


// spatial filter - a function we'll use inside our synthDef below...
spatFilt = { arg in, coeffMatrix;

	// wrap input as array if needed, for mono inputs
	in.isArray.not.if({ in = [in] });

	Mix.fill( coeffMatrix.cols, { arg i; // fill input
		UGen.replaceZeroesWithSilence(
			coeffMatrix.asArray.reshape(coeffMatrix.rows, coeffMatrix.cols).flop.at(i) * in.at(i)
		)
	})
};

// a function to return coefficients for Balance- Matrix form
// position argument in degrees
balanceMatrix = { arg angleInDegrees = 0;
	var angleInRadians;
	var theta;

	angleInRadians = angleInDegrees/180*pi;

	theta = pi/4 - angleInRadians;

	Array2D.fromArray(
		2, // rows, outputs
		2, // columns, inputs
		2.sqrt * [
			theta.cos, 0,
			0, theta.sin
		]
	)
};

// a function to return coefficients for Width- Matrix form
// width argument in degrees
widthMatrix = { arg angleInDegrees = 0;
	var angleInRadians;
	var theta;

	angleInRadians = angleInDegrees/180*pi;

	theta = angleInRadians;

	Array2D.fromArray(
		2, // rows, outputs
		2, // columns, inputs
		[
			theta.cos, -1 * theta.sin,
			-1 * theta.sin, theta.cos
		]
	)
};

///////////////    DEFINE FILTER FUNCTIONS    ///////////////

sinc = { arg x;
	(x == 0).if({ 1 }, { sin(pi * x) / (pi * x) })
};

lopass = { arg size, freq, sampleRate;
	var bandWidth;

	// compute normalised bandwidth
	bandWidth = 2 * freq / sampleRate;

	// generate kernel
	size.collect({arg n;
		bandWidth * sinc.value(bandWidth * (n - ((size-1)/2)))
	}).as(Signal)
};

hipass = { arg size, freq, sampleRate;
	lopass.value(size, sampleRate/2, sampleRate) - lopass.value(size: size, freq: freq, sampleRate: sampleRate);
};

// sinc bandpass - undefined!
bapass = { arg size, freq, bw, sampleRate;
	var freq0, freq1;

	// compute freqs
	freq0 = freq - (bw/2);
	freq1 = freq + (bw/2);

	lopass.value(size, freq1, sampleRate) - lopass.value(size, freq0, sampleRate); // LP1 - LP0
};

// sinc bandstop (complementary)
bastop = { arg size, freq, bw, sampleRate;
	lopass.value(size, sampleRate/2, sampleRate) - bapass.value(size, freq, bw, sampleRate); // AP - BP
};

// sinc low-shelf (complementary)
loshelf = { arg size, freq, gain, sampleRate;
	var scale = gain.dbamp;
	hipass.value(size, freq, sampleRate) + (lopass.value(size, freq, sampleRate) * scale);
};

// sinc high-shelf (complementary)
hishelf = { arg size, freq, gain, sampleRate;
	var scale = gain.dbamp;
	(hipass.value(size, freq, sampleRate) * scale) + lopass.value(size, freq, sampleRate);
};


// damping as a function of order, n, and stage, k
dampFacFunc = { arg n, k;
	var d;

	d = sin((((2*(k + 1)) - 1) * pi) / (2 * n));

	d;
};

// HP || LP filter coefficients - SOS
sosFiltCoeffsFunc = { arg func, freq, d, sr;
	var c;
	var a0, a1, a2, b1, b2;

	switch(func,
		1, {
			c = tan(pi * freq / sr);

			a0 = 1/(1 + (2 * d * c) + (c.squared));
			a1 = -2 * a0;
			a2 = a0;
			b1 = 2 * a0 * (1 - (c.squared));
			b2 = -1 * a0 * (1 - (2 * d * c) + (c.squared));
		},
		2, {
			c = 1/tan(pi * freq / sr);

			a0 = 1/(1 + (2 * d * c) + (c.squared));
			a1 = 2 * a0;
			a2 = a0;
			b1 = -1 * 2 * a0 * (1 - (c.squared));
			b2 = -1 * a0 * (1 - (2 * d * c) + (c.squared));
		}
	);

	Array.with(a0, a1, a2, b1, b2);
};


// calculates coeffs for HP (func = 1) || LP (func = 2) filters of order n
cascadeCoeffsFunc = { arg func, freq, n, sr;
	((n/2).asInt).collect({ arg k;
		sosFiltCoeffsFunc.value(
			func,
			freq,
			dampFacFunc.value(n, k),
			sr
		)
	})
};

// Schroeder Allpass (Even)
schroederCoeffFunc = { arg bFac;

	var tanFac, gFac;

	tanFac = tan(pi/2 * bFac);
	gFac = (1 - tanFac) / (1 + tanFac);

	// return
	gFac;
};

// function to calculate Regalia-Mitra shelf coefficient
//
// Regalia-Mitra
regaliaMitraCoeff = { arg gain;

	var kFac;

	kFac = (gain.dbamp - 1) / 2;

	// return
	kFac;
};

// function to calculate delay
//
// Comb Filter Delay
combDelayFunc = { arg freq;

	var delay;

	delay = (2*freq).reciprocal;

	// return
	delay;
};

/////////////// DEFINE SYNTH HELPER FUNCTIONS ///////////////

// sine-cosine panning law coefficient function
// angle argument in degrees
sinCosPanLaw = { arg angleInDegrees = 0;
	var angleInRadians;
	var theta;

	angleInRadians = angleInDegrees/180*pi;

	theta = pi/4 - angleInRadians;

	[theta.cos, theta.sin]
};

equalPowMix= { arg mix = 0.5;

	[(1 - mix).sqrt, mix.sqrt]
};

// converts from Dattoro Sample Rate of 29761Hz to this sampleRate
sampConvert = {arg num;
	(num / 29761) * sampleRate
};

allPoleCoeffFunc = { arg delayTime, decayTime;

	var gFac;

	gFac = 10.pow(-3 * delayTime / decayTime);

	// return
	gFac;
};

combT60Func = { arg delay, gFac;

	var t60;

	t60 = gFac.sign * (-3 * delay / log10(gFac.abs));

	// return
	t60;
};

sampsToTime = { arg numSamples, sampleRate;

	var time;

	time = numSamples / sampleRate;

	// return
	time;
};

// -1 = first channel, 1 = second channel
equalPowMix = { arg mix = 0.0;
	mix = mix + 1;
	mix = mix / 2;
	[(1 - mix).sqrt, mix.sqrt]
};

// Sinus Delay Time
vibSinDelTimeFunc = { arg ratio, rate;

	var delayTime;

	delayTime = (ratio - 1) / (pi * rate);

	// return
	delayTime;
};

///////////////// DEFINE SYNTHDEFS //////////////////

/* ------------- PLAYBACK ------------- */

// basic playback synth
// dur = duration of buffer, gain = gain applied to output, ampEnv = amp envelope, rate = playback rate, loop = loop flag,
// buffer = buffer for playback, sendBus1 = output bus, sendBus2 = output bus
// startPos = start position of playback
playBufSynthDef = CtkSynthDef.new(\myStereoPlayBufSynth, {arg dur, gain = 0.0, startPos = 0, ampEnv = 1.0, rate = 1, loop = 0, buffer = 0, sendBus1 = 0.0, sendBus2 = 0.0;

	var numChannels = 2; // stereo buffer
	var amp, sendAmp;          // a few vars for synthesis
	var sig;     // vars assigned to audio signals
	var out;

	// calcs
	amp = gain.dbamp;  // convert from gain in dB to linear amplitude scale

	// read-in signal & envelope
	sig = PlayBuf.ar(numChannels, buffer,  BufRateScale.kr(buffer) * rate, loop: loop, startPos: startPos * SampleRate.ir);

	// apply amp envelope
	sig = ampEnv * sig;

	out = amp * sig;

	Out.ar(sendBus1, out);

	Out.ar(sendBus2, out);
});


// Synchronous Granular Synthesis for independent time stretching and pitch shifting
sgsSynth = CtkSynthDef.new(\sgsSynth, {arg dur, gain, bufL = 0, bufR = 0, startPos = 0.0, endPos = 1.0, ampEnv = 1, freq = 440, wavRatio = 1.0, refFreq = 440, q = 1.0, sendBus1 = 0, sendBus2 = 0;

	var numChannels = 1; // num channels read by GrainBuf
	var grainFreq, envFreq, grainDur, bufferPersL, bufferPersR;
	var amp;
	var indx, indxL, indxR, trigger, sig, sigL, sigR;

	// calcs
	amp = gain.dbamp;
	grainFreq = freq;
	envFreq = wavRatio * refFreq / (2 * q);
	grainDur = envFreq.reciprocal;
	bufferPersL = refFreq * BufDur.kr(bufL);
	bufferPersR = refFreq * BufDur.kr(bufR);

	// buffer position (index)
	indx = Line.ar(startPos, endPos, dur);
	indxL = (indx * bufferPersL).floor / bufferPersL;
	indxR = (indx * bufferPersR).floor / bufferPersR;

	trigger = Impulse.ar(grainFreq);

	// granular synthesis & envelope scaling
	// insert two channel processing
	sigL = GrainBuf.ar(numChannels: numChannels, trigger: trigger, dur: grainDur, sndbuf: bufL, rate: wavRatio, pos: indxL);
	sigR = GrainBuf.ar(numChannels: numChannels, trigger: trigger, dur: grainDur, sndbuf: bufR, rate: wavRatio, pos: indxR);
	sig = [sigL, sigR];
	sig = ampEnv * amp * sig;

	// out!!
	Out.ar(sendBus1, sig);

	Out.ar(sendBus2, sig);

});

// Asynchronous Granular Synthesis for independent time stretching and pitch shifting
agsSynth = CtkSynthDef.new(\agsSynth, {arg dur, gain, bufL = 0, bufR = 0, startPos = 0.0, endPos = 1.0, ampEnv = 1, density = 1, wavRatio = 1.0, refFreq = 440, q = 1.0, sendBus1 = 0, sendBus2 = 0;

	var numChannels = 1; // num channels read by GrainBuf
	var grainFreq, envFreq, grainDur, bufferPersL, bufferPersR;
	var amp;
	var indx, indxL, indxR, trigger, sig, sigL, sigR;

	// calcs
	amp = gain.dbamp;
	grainFreq = density * q.reciprocal * wavRatio * refFreq;
	envFreq = wavRatio * refFreq / (2 * q);
	grainDur = envFreq.reciprocal;
	bufferPersL = refFreq * BufDur.kr(bufL);
	bufferPersR = refFreq * BufDur.kr(bufR);

	// buffer position (index)
	indx = Line.ar(startPos, endPos, dur);
	indxL = (indx * bufferPersL).floor / bufferPersL;
	indxR = (indx * bufferPersR).floor / bufferPersR;

	trigger = Dust.ar(grainFreq);

	// granular synthesis & envelope scaling
	// insert two channel processing
	sigL = GrainBuf.ar(numChannels: numChannels, trigger: trigger, dur: grainDur, sndbuf: bufL, rate: wavRatio, pos: indxL);
	sigR = GrainBuf.ar(numChannels: numChannels, trigger: trigger, dur: grainDur, sndbuf: bufR, rate: wavRatio, pos: indxR);
	sig = [sigL, sigR];
	sig = ampEnv * amp * sig;

	// out!!
	Out.ar(sendBus1, sig);

	Out.ar(sendBus2, sig);

});

/* ------------- FILTERS ------------- */

// hp = highpass, lp = lowpass, bp = bandpass, br = bandreject
sosHPFilter = CtkSynthDef.new(\sosHPFilter, {arg freq = 440, receiveBus, sendBus = 0;
	var sig, coeffs;
	var numChannels = 2;
	var order = 32;

	sig = In.ar(receiveBus, numChannels);

	coeffs = cascadeCoeffsFunc.value(1, freq, order, SampleRate.ir);

	coeffs.do({arg stage;
		sig = SOS.ar(sig, stage.at(0), stage.at(1), stage.at(2), stage.at(3), stage.at(4));
	});

	Out.ar(sendBus, sig);
});

sosLPFilter = CtkSynthDef.new(\sosLPFilter, {arg freq = 440, receiveBus, sendBus = 0;
	var sig, coeffs;
	var numChannels = 2;
	var order = 32;

	sig = In.ar(receiveBus, numChannels);

	coeffs = cascadeCoeffsFunc.value(2, freq, order, SampleRate.ir);

	coeffs.do({arg stage;
		sig = SOS.ar(sig, stage.at(0), stage.at(1), stage.at(2), stage.at(3), stage.at(4));
	});

	Out.ar(sendBus, sig);

});

sosBPFilter = CtkSynthDef.new(\sosBPFilter, {arg freq = 440, bWidth = 10, receiveBus, sendBus = 0;
	var sig, coeffs1, coeffs2;
	var numChannels = 2;
	var order = 32;

	sig = In.ar(receiveBus, numChannels);

	coeffs1 = cascadeCoeffsFunc.value(1, freq - (bWidth / 2), order, SampleRate.ir);

	coeffs2 = cascadeCoeffsFunc.value(2, freq + (bWidth / 2), order, SampleRate.ir);

	coeffs1.do({arg stage;
		sig = SOS.ar(sig, stage.at(0), stage.at(1), stage.at(2), stage.at(3), stage.at(4));
	});

	coeffs2.do({arg stage;
		sig = SOS.ar(sig, stage.at(0), stage.at(1), stage.at(2), stage.at(3), stage.at(4));
	});

	Out.ar(sendBus, sig);
});

sosBSFilter = CtkSynthDef.new(\sosBSFilter, {arg freq = 440, bWidth = 10, receiveBus, sendBus = 0;
	var sig, coeffs1, coeffs2;
	var numChannels = 2;
	var order = 32;

	sig = In.ar(receiveBus, numChannels);

	coeffs1 = cascadeCoeffsFunc.value(2, freq - (bWidth / 2), order, SampleRate.ir);

	coeffs2 = cascadeCoeffsFunc.value(1, freq + (bWidth / 2), order, SampleRate.ir);

	coeffs1.do({arg stage;
		sig = SOS.ar(sig, stage.at(0), stage.at(1), stage.at(2), stage.at(3), stage.at(4));
	});

	coeffs2.do({arg stage;
		sig = SOS.ar(sig, stage.at(0), stage.at(1), stage.at(2), stage.at(3), stage.at(4));
	});

	Out.ar(sendBus, sig);
});

firFilter = CtkSynthDef.new(\firFilter, {arg buffer = 0.0, receiveBus, sendBus = 0;

	var sig;
	var numChannels = 2;

	sig = In.ar(receiveBus, numChannels);

	sig = Convolution2.ar(sig, buffer, framesize: buffer.size);

	Out.ar(sendBus, sig);

});

convBasic = CtkSynthDef.new(\convBasic, {arg buffer, receiveBus, sendBus = 0;
	var sig;
	var numChannels = 2;

	sig = In.ar(receiveBus, numChannels);
	sig = Convolution2.ar(sig, buffer, framesize: buffer.size);

	Out.ar(sendBus, sig);
});


// unimplemented, intended to use one buffer as a spectral envelope for another
// possibly use (ctk note instance).createBuffer -> buffer frame to kernel
convEnv = CtkSynthDef.new(\convBasic, {arg receiveBusA, bufferA, dur, bufferSamp, sendBus = 0;
	/* var sigA, sigB, sig;
	var numChannels = 1;
	var indx, kernel;

	sigA = In.ar(receiveBusA, numChannels);

	indx = Line.ar(0.001, 1, dur);

	kernel =

	sig = PartConv.ar(sigA, 1024, kernel);

	Out.ar(sendBus, sig); */
});


// implements even || odd comb filter, even = 1 || even = 0
// freq = frequency to pass at even or odd harmonics
// minFreq = minimum freq to pass
// gain = how much to reject or boost non-passed freqs
// b = bandwidth of comb
universalComb = CtkSynthDef.new(\evenComb, { arg even = 1, freq = 440.0, minFreq = 20.0, combGain = -60.0, b = 0.5, receiveBus, sendBus = 0;
	var sig;
	var impulse;
	var maxDelayTime, delayTime, decayTime;
	var gFac, kFac;
	var numChannels = 2;

	// comb filter calcs
	maxDelayTime = combDelayFunc.value(minFreq);
	delayTime = combDelayFunc.value(freq);
	gFac = schroederCoeffFunc.value(b);
	decayTime = combT60Func.value(delayTime, gFac);
	kFac = regaliaMitraCoeff.value(combGain);

	// generate impulse
	impulse = In.ar(receiveBus, numChannels);

	// test filter
	sig = impulse + (kFac * AllpassC.ar(impulse, maxDelayTime, delayTime, (-1 + (2 * even)) * decayTime, (-2 * even) + 1, impulse));

	Out.ar(sendBus, sig);
});

/* ------------- SPECIAL FX ------------- */

// generates feedback echo only (doesn't include source)
feedbackEcho = CtkSynthDef.new(\feedbackSynth, {arg gain = 0.0, ampEnv = 1.0, delayTime = 0.2, maxDelayTime = 0.2, decayTime = 0.5, receiveBus, sendBus = 0;
	var numChannels = 2;
	var amp;
	var sig, delay, out;     // vars assigned to audio signals

	// calcs
	amp = gain.dbamp;

	// read sound in
	sig = In.ar(receiveBus, numChannels);

	// delay
	delay = CombN.ar(sig, maxDelayTime, delayTime, decayTime);

	// apply gain to echo effect
	out = delay * amp * ampEnv;

	// out!!
	Out.ar(sendBus, out)
});

// generates feedback echo only (doesn't include source)
slapbackEcho = CtkSynthDef.new(\slapbackSynth, {arg gain = 0.0, ampEnv = 1.0, delayTime = 0.2, maxDelayTime = 0.2, receiveBus, sendBus = 0;
	var numChannels = 2;
	var amp;
	var sig, delay, out;     // vars assigned to audio signals

	// calcs
	amp = gain.dbamp;

	// read sound in
	sig = In.ar(receiveBus, numChannels);

	// delay
	delay = DelayN.ar(sig, maxDelayTime, delayTime);

	// apply gain to echo effect
	out = delay * amp * ampEnv;

	// out!!
	Out.ar(sendBus, out)
});


// Dattorro Tank Reverb
// outputs source -> mix = -1.0, reverb mix = 1.0
// modified to also accept a custom filter kernel for versatility
dReverb = CtkSynthDef.new(\dReverb, {arg dur, gain = 0.0, ampEnv = 1, preDelayTime = 0.2, maxDelayTime = 0.2, lowCut = 200.0, filter = 1, receiveBus, sendBus = 0;
	var numChannels = 2; // stereo bus!
	var amp, dirAmp, revAmp;          // a few vars for synthesis
	var direct, out;     // vars assigned to audio signals
	var del, time; // temp vars

	// variables for local/private bussing
	var sigA, sigB;
	var reverb;

	// variables for output taps
	var yL = 0;
	var yR = 0;
	var l48, l54, l59, l63, l30, l33, l39;
	var r24, r30, r33, r39, r54, r59, r63;

	// variables for reverberation
	var sampR = 29761;
	var excurs = 8; // half of original value for simpler calcs using LFNoise2
	var revDecay = 0.5;
	var bandwidth = 0.9995; // QUESTION: what is this parameter's purpose?
	var damping = 0.0005; // QUESTION: what is this parameter's purpose?
	var sapCoeffs = [142, 107, 379, 277];
	var inputDiff = [0.75, 0.75, 0.5, 0.5];
	var sigAParams = [[(908 + excurs) + (LFNoise2.kr(1) * excurs), 0.70, 4217, 1 - damping, lowCut, revDecay], [2656, 0.50, 3163, 1, sampleRate / 2, revDecay]]; // params from Dattarro ugen graph, some excluded
	var sigBParams = [[(672 + excurs) + (LFNoise2.kr(1) * excurs), 0.70, 4453, 1 - damping, lowCut, revDecay], [1800, 0.50, 3720, 1, sampleRate / 2, revDecay]]; // params from Dattarro ugen graph, some excluded

	// calcs
	amp = gain.dbamp;  // overall scalar

	// read sound in
	direct = In.ar(receiveBus, numChannels);

	// multiply signal
	reverb = direct * 0.5;

	// apply delay (delayTime)
	reverb = DelayC.ar(in: reverb, maxdelaytime: maxDelayTime, delaytime: preDelayTime, mul: 1, add: 0);

	// apply lowpass filter (lowCut)
	reverb = LPF.ar(in: reverb, freq: lowCut, mul: 1, add: 0);
	reverb = Convolution2.ar(reverb, filter, framesize: filter.size);

	// loop for SAP1, SAP2, (revDecay & input diff1 = 0.750) SAP3, SAP4 (revDecay & input diff2 = 0.625)
	sapCoeffs.do({arg c, i;
		var d = sampsToTime.value(sapCoeffs.at(i), sampR);
		var t = combT60Func.value(d, inputDiff.at(i));
		reverb = AllpassC.ar(reverb, maxDelayTime, d, t);
	});

	// SIG A: add SigB, SAP5 (w/ excursion), delay, LPF, decay multiplier, SAP6, delay, decay multiplier
	sigA = reverb + LocalIn.ar(numChannels)[0];

	del = sampsToTime.value(sigAParams.at(0).at(0), sampR);
	time = combT60Func.value(del, sigAParams.at(0).at(1));
	sigA = AllpassC.ar(sigA, maxDelayTime, del, time); // all pass
	l48 = sigA; // LEFT TAP
	sigA = DelayC.ar(sigA, maxDelayTime, sampsToTime.value(sigAParams.at(0).at(2), sampR)); // delay
	l54 = sigA; // LEFT TAP
	r54 = sigA; // RIGHT TAP
	sigA = LPF.ar(in: sigA, freq: sigAParams.at(0).at(4), mul: 1, add: 0);
	sigA = Convolution2.ar(sigA, filter, framesize: filter.size);
	sigA = sigA * sigAParams.at(0).at(5);

	del = sampsToTime.value(sigAParams.at(1).at(0), sampR);
	time = combT60Func.value(del, sigAParams.at(1).at(1));
	sigA = AllpassC.ar(sigA, maxDelayTime, del, time); // all pass
	l59 = sigA; // LEFT TAP
	r59 = sigA; // RIGHT TAP
	sigA = DelayC.ar(sigA, maxDelayTime, sampsToTime.value(sigAParams.at(1).at(2), sampR)); // delay
	l63 = sigA; // LEFT TAP
	r63 = sigA; // RIGHT TAP
	sigA = sigA * sigAParams.at(1).at(5);

	// SIG B: add SigA, SAP7 (w/ excursion), delay, LPF, decay multiplier, SAP8, delay, decay multiplier
	sigB = reverb + LocalIn.ar(numChannels)[1];

	del = sampsToTime.value(sigBParams.at(0).at(0), sampR);
	time = combT60Func.value(del, sigBParams.at(0).at(1));
	sigB = AllpassC.ar(sigB, maxDelayTime, del, time); // all pass
	r24 = sigB; // RIGHT TAP
	sigB = DelayC.ar(sigB, maxDelayTime, sampsToTime.value(sigBParams.at(0).at(2), sampR)); // delay
	l30 = sigB; // LEFT TAP
	r30 = sigB; // RIGHT TAP
	sigB = LPF.ar(in: sigB, freq: sigBParams.at(0).at(4), mul: 1, add: 0);
	sigB = Convolution2.ar(sigB, filter, framesize: filter.size);
	sigB = sigB * sigBParams.at(0).at(5);

	del = sampsToTime.value(sigBParams.at(1).at(0), sampR);
	time = combT60Func.value(del, sigBParams.at(1).at(1));
	sigB = AllpassC.ar(sigB, maxDelayTime, del, time); // all pass
	l33 = sigB; // LEFT TAP
	r33 = sigB; // RIGHT TAP
	sigB = DelayC.ar(sigB, maxDelayTime, sampsToTime.value(sigBParams.at(1).at(2), sampR)); // delay
	l39 = sigB; // LEFT TAP
	r39 = sigB; // RIGHT TAP
	sigB = sigB * sigAParams.at(1).at(5);

	// Feedback sigA and sigB
	LocalOut.ar([sigA, sigB]);

	// add taps for output
	yL = (0.6)*(l48  + l54 - l59 + l63 - l30 - l33 - l39);
	yR = (0.6)*(r39  + r30 - r24 + r39 - r54 - r59 - r63);

	// add panned direct and delay to out & envelope
	out = amp * ampEnv * [yL, yR];

	// out!!
	Out.ar(sendBus, out)
});

// vibrato & tremolo using LFNoise to modulate delay times
vibTrem = CtkSynthDef.new(\vibTrem, {arg gain = 0, ampEnv = 1, rate = 6.0, minRate = 6.0, ratio = 1.01, maxRatio = 1.01, modIndexGain = -12.0, changeRate = 10, change = 0.05, receiveBus, sendBus = 0;

	var amp, sig, delay, out;
	var maxDelayTime, delayTime, modIndex, normFac, modLFO;
	var numChannels = 2;

	// calcs
	amp = gain.dbamp;  // convert from gain in dB to linear amplitude scale
	maxDelayTime = vibSinDelTimeFunc.value(maxRatio, minRate);
	delayTime = vibSinDelTimeFunc.value(ratio, rate);
	maxDelayTime = maxDelayTime + (maxDelayTime * LFNoise2.kr(changeRate, change));
	delayTime = delayTime + (delayTime * LFNoise2.kr(changeRate, change));

	modIndex = modIndexGain.dbamp; // convert AM gain in dB to linear amplitude scale
	normFac = (1 + (2*modIndex)).reciprocal; // AM amplitude normalization factor

	sig = In.ar(receiveBus, numChannels);

	// delay line & amplitude modulator (unscaled)
	modLFO = SinOsc.ar(rate, [pi/2, 0]); // quadrature modulator

	// modulated delay line
	delay = DelayC.ar(sig, maxDelayTime, modLFO.at(0).range(0, 1) * delayTime);  // vibrato
	delay = normFac * (1 + (modIndex * modLFO.at(1))) * delay; // tremelo

	// add panned direct and delay to out & envelope
	out = amp * ampEnv * delay;

	// out!!
	Out.ar(sendBus, out);
});


// weaverSSB in STEREO
// ssbFreq
// bw
weaverSSB = CtkSynthDef.new(\weaverSSB, { arg gain = 0.0, ampEnv = 1, ssbFreq = 1.0, bw = 1.0, receiveBus, sendBus = 0;

	// variables
	var numChannels = 2; // stereo
	var srDiv4, quadOsc0, quadOsc1;     // vars assigned to audio signals
	var sig, sigL, sigR, ssbSigL, ssbSigR, ssbSig, out;
	var amp;          // a few vars for synthesis
	var cascadeCoeffs;
	// var n = 2; // filter order
	// var n = 4; // filter order
	// var n = 6; // filter order
	// var n = 8; // filter order
	// var n = 16; // filter order
	// var n = 24; // filter order
	// var n = 32; // filter order
	var n = 64; // filter order


	// calcs
	amp = gain.dbamp;
	srDiv4 = SampleRate.ir/4; // SSB calcs
	cascadeCoeffs = cascadeCoeffsFunc.value(func: 2, freq: bw * srDiv4, n: n, sr: SampleRate.ir); // low pass filter calcs

	// the sample playback oscillator
	sig = In.ar(receiveBus, numChannels);

	// seperate left and right channels
	sigL = sig[0];

	sigR = sig[1];

	// Weaver SSB - quadrature oscillators
	quadOsc0 = SinOsc.ar(srDiv4, [pi/2, 0]);
	quadOsc1 = SinOsc.ar(srDiv4 + ssbFreq, [pi/2, 0]);

	// Weaver SSB - lowpass filter
	ssbSigL = sigL * quadOsc0;  // <-- input to the Weaver Modulation Network here!
	ssbSigR = sigR * quadOsc0;
	cascadeCoeffs.do({ arg coeffs;
		ssbSigL = SOS.ar(ssbSigL, coeffs.at(0), coeffs.at(1), coeffs.at(2), coeffs.at(3), coeffs.at(4));
		ssbSigR = SOS.ar(ssbSigR, coeffs.at(0), coeffs.at(1), coeffs.at(2), coeffs.at(3), coeffs.at(4));
	});
	ssbSigL = ssbSigL * quadOsc1;
	ssbSigR = ssbSigR * quadOsc1;
	ssbSig = [ssbSigL, ssbSigR];
	ssbSig = 2 * ssbSig;  // restore gain - NOTE: this could be applied elsewhere

	// envelope
	ssbSig = amp * ampEnv * ssbSig;

	// out!!
	Out.ar(sendBus, ssbSig)
});

// absolute value shaper
absShaper = CtkSynthDef.new(\absShaper, {arg receiveBus, sendBus = 0;
	var sig;
	var numChannels = 2;

	sig = In.ar(receiveBus, numChannels);

	Out.ar(sendBus, sig.abs);
});


// cos shaper
cosShaper = CtkSynthDef.new(\cosShaper, {arg receiveBus, sendBus = 0;
	var sig;
	var numChannels = 2;

	sig = In.ar(receiveBus, numChannels);

	Out.ar(sendBus, sig.cos);
});

// distort
disShaper = CtkSynthDef.new(\cosShaper, {arg receiveBus, sendBus = 0;
	var sig;
	var numChannels = 2;

	sig = In.ar(receiveBus, numChannels);

	Out.ar(sendBus, sig.distort);
});

//

// add in custom distortion functions (sin, cos, tan, abs etc)
// add in SSB here
// try to emulate distortion from lecture code using bussing principles

/* ------------- IMAGING ------------- */

//

// designed such that the client has the choice to keep busses seperate but simulate crossfade
crossFade = CtkSynthDef.new(\crossFader, {arg fade = 0.0, receiveBusA, receiveBusB, sendBusA = 0, sendBusB = 0;

	var sigA, sigB, out;
	var numChannels = 2;

	sigA = In.ar(receiveBusA, numChannels);
	sigB = In.ar(receiveBusB, numChannels);
	sigA = sigA * equalPowMix.value(fade).at(0);
	sigB = sigB * equalPowMix.value(fade).at(1);
	Out.ar(sendBusA, sigA);
	Out.ar(sendBusB, sigB);
});

rotateImage = CtkSynthDef.new(\rotateImage, {arg angle = 0, receiveBus, sendBus = 0;
	var sig;
	var numChannels = 2;

	sig = In.ar(receiveBus, numChannels);

	// spatial filter
	sig = spatFilt.value(sig, rotateMatrix.value(angle));  // <-- stereo imaging happens here!

	// out!!
	Out.ar(sendBus, sig)
});

balanceImage = CtkSynthDef.new(\balanceImage, {arg angle = 0, receiveBus, sendBus = 0;
	var sig;
	var numChannels = 2;

	sig = In.ar(receiveBus, numChannels);

	// spatial filter
	sig = spatFilt.value(sig, balanceMatrix.value(angle));  // <-- stereo imaging happens here!

	// out!!
	Out.ar(sendBus, sig)
});

widthImage = CtkSynthDef.new(\widthImage, {arg angle = 0, receiveBus, sendBus = 0;
	var sig;
	var numChannels = 2;

	sig = In.ar(receiveBus, numChannels);

	// spatial filter
	sig = spatFilt.value(sig, widthMatrix.value(angle));  // <-- stereo imaging happens here!

	// out!!
	Out.ar(sendBus, sig)
});


// comb filter frequency spreader implemented w/ universal comb, stereo in -> stereo out
// evenLeft = 1 -> left channel has even combs, right channel has odd combs
// evenLeft = 0 -> left channel has odd combs, right channel has even combs
// freq = frequency to pass at even or odd harmonics
// minFreq = minimum freq to pass
// gain = how much to reject or boost non-passed freqs
// b = bandwidth of comb
combSpreader = CtkSynthDef.new(\combSpreader, {arg evenLeft = 1, freq = 440.0, minFreq = 20.0, combGain = -60, b = 0.5,  widthAngle = 0.0, receiveBus, sendBus = 0;

	// variables
	var numChannels = 2; // stereo buffer
	var sig, sigL, sigR, combL, combR, out;     // vars assigned to audio signals
	var maxDelayTime, delayTime;
	var gFac, decayTime, kFac;
	var evenRight = ((-1 + evenLeft) * -1); // evenLeft = 1 -> evenRight = 0, evenLeft = 0 -> evenRight = 1

	// comb filter calcs
	maxDelayTime = combDelayFunc.value(minFreq);
	delayTime = combDelayFunc.value(freq);
	gFac = schroederCoeffFunc.value(b);
	decayTime = combT60Func.value(delayTime, gFac);
	kFac = regaliaMitraCoeff.value(combGain);

	// the sample playback oscillator - mono
	sig = In.ar(receiveBus, numChannels);

	// complementary combs - cos & sin
	combL = sig + (kFac * AllpassC.ar(sig, maxDelayTime, delayTime, (-1 + (2 * evenLeft)) * decayTime, (-2 * evenLeft) + 1, sig));  // comb for left channel
	combR = sig + (kFac * AllpassC.ar(sig, maxDelayTime, delayTime, (-1 + (2 * evenRight)) * decayTime, (-2 * evenRight) + 1, sig)); // comb for right channel

	// spatial filter - width & rotate
	combL = spatFilt.value(combL, rotateMatrix.value(widthAngle));  // creates width spreading
	combR = spatFilt.value(combR, rotateMatrix.value(-1 * widthAngle));

	// envelope scaling
	out = combL + combR;

	// out!!
	Out.ar(sendBus, out)
});


// simulates speaker rotating at variable rate, sending left and right channels out different cones at
// a specific angle
rotarySpreader = CtkSynthDef.new(\rotarySpreader, {arg gain = 0, ampEnv = 1, rate = 6.0, minRate = 6.0, ratio = 1.01, maxRatio = 1.01, modIndexGain = -12.0, changeRate = 10, change = 0.05, rotateAngle = 0.0, receiveBus, sendBus = 0;

	var numChannels = 2; // mono bus!
	var maxDelayTime, delayTime;
	var dry, sigL, sigR, out, amp;     // vars assigned to audio signals
	var delay;
	var modLFO;
	var normFac;      // normalization factor
	var modIndex;     // modulation index (a scalar)

	// calcs
	amp = gain.dbamp;  // convert from gain in dB to linear amplitude scale
	maxDelayTime = vibSinDelTimeFunc.value(maxRatio, minRate);
	delayTime = vibSinDelTimeFunc.value(ratio, rate);
	maxDelayTime = maxDelayTime + (maxDelayTime * LFNoise2.kr(changeRate, change));
	delayTime = delayTime + (delayTime * LFNoise2.kr(changeRate, change));

	modIndex = modIndexGain.dbamp; // convert AM gain in dB to linear amplitude scale
	normFac = (1 + (2*modIndex)).reciprocal; // AM amplitude normalization factor

	// read sound in
	dry = In.ar(receiveBus, numChannels);

	sigL = dry[0]; // left channel

	sigR = dry[1]; // right channel

	// LEFT CHANNEL MODULATION
	// delay line & amplitude modulator (unscaled)
	modLFO = SinOsc.ar(rate, [pi/2, 0]); // quadrature modulator, normal phase

	// modulated delay line
	sigL = DelayC.ar(sigL, maxDelayTime, modLFO.at(0).range(0, 1) * delayTime);  // vibrato
	sigL = normFac * (1 - (modIndex * modLFO.at(1))) * sigL; // tremelo


	// RIGHT CHANNEL MODULATION

	// modulated delay line
	sigR = DelayC.ar(sigR, maxDelayTime, (-1 * modLFO.at(0).range(0, 1)) * delayTime);  // vibrato
	sigR = normFac * (1 + (modIndex * modLFO.at(1))) * sigR; // tremelo

	// add panned direct and delay to out & envelope
	out = [sigL, sigR]; // add chanels together
	out = spatFilt.value(out, rotateMatrix.value(rotateAngle)); // unequal mixing using rotation
	out = amp * ampEnv * out;

	// out!!
	Out.ar(sendBus, out)
});

/////////////// DEFINE NOTE CALLING FUNCTIONS ///////////////

getLeft = {arg buffer;

	var left = CtkBuffer.new(buffer.path, numChannels: 1, channels: 0);
	score.add(left);
	left;
};

getRight = {arg buffer;

	var right = CtkBuffer.new(buffer.path, numChannels: 1, channels: 1);
	score.add(right);
	right;

};


///////////////// CREATE BUFFERS //////////////////

// get path of sampled sounds folder
sndPath = "continious/".resolveRelative;
// get individual file paths
allSounds = (sndPath ++ "*").pathMatch;

// collects Ctk buffers into a dictionary with keys of filenames
buffersC = Dictionary.new;
allSounds.do{arg bufferPath;
	var name;
	bufferPath = bufferPath.replace("\\", "/"); // converts windows path to unix path
	name = bufferPath.basename.splitext.at(0); // gets just the file name
	buffersC.add((name) -> (CtkBuffer.new(bufferPath)));
	score.add(buffersC.at(name));
};

// get path of sampled sounds folder
sndPath = "discrete/".resolveRelative;
// get individual file paths
allSounds = (sndPath ++ "*").pathMatch;

// collects Ctk buffers into a dictionary with keys of filenames
buffersD = Dictionary.new;
allSounds.do{arg bufferPath;
	var name;
	bufferPath = bufferPath.replace("\\", "/"); // converts windows path to unix path
	name = bufferPath.basename.splitext.at(0); // gets just the file name
	buffersD.add((name) -> (CtkBuffer.new(bufferPath)));
	score.add(buffersD.at(name));
};

///////////////// POPULATE THE SCORE //////////////////

// add the buffer and node group to the score
// NOTE: buffers & node groups must be added to the score for the CtkSynthDef to access!
score.add(group);

// basic envelopes and filters
allPassBuffer = CtkBuffer.collection(lopass.value(1024, sampleRate / 2, sampleRate));
score.add(allPassBuffer);

basicAmpEnv = Env([0, 1, 1, 0], [0.05, 0.9, 0.05], \lin);

// #### INTRO VOCALS

buffer = buffersD.at("vocal_talking_motion_1");
start = 0.0;
dur = buffer.duration;
score.add(
	sgsSynth.note(starttime: start, duration: dur, addAction: \head, target: group)
	.dur_(dur)
	.wavRatio_(CtkControl.env(Env([0.5, 1.0], [2.3], \lin), starttime: start))
	.refFreq_(100)
	.freq_(CtkControl.env(Env([50, 100], [2.3], \lin), starttime: start))
	.gain_(12.0)
	.bufL_(getLeft.value(buffer))
	.bufR_(getRight.value(buffer))
	.sendBus1_(channels[0][0])
	.sendBus2_(garbage)
);

// reverb
score.add(
	dReverb.note(starttime: start + 4.0, duration: 2.0, addAction: \tail, target: group)
	.filter_(allPassBuffer)
	.receiveBus_(channels[0][0])
	.sendBus_(0)
);

// width
score.add(
	widthImage.note(starttime: start, duration: dur, addAction: \tail, target: group)
	.angle_(CtkControl.env(Env([0, 0, 45, 45, 0], [7.11, 0.01, 2.3, 0.01], \step), starttime: start))
	.receiveBus_(channels[0][0])
	.sendBus_(0)
);

// vibrato/tremolo
score.add(
	vibTrem.note(starttime: start + 13, duration: 5.0, addAction: \tail, target: group)
	.modIndexGain_(CtkControl.env(Env([0, 0, -15.0], [13, 0.01], \step), starttime: start))
	.rate_(4.0)
	.minRate_(4.0)
	.ratio_(1.02)
	.maxRatio_(1.05)
	.receiveBus_(channels[0][0])
	.sendBus_(0)
);

// echo out
score.add(
	feedbackEcho.note(starttime: start + 15.4, duration: 3.5, addAction: \tail, target: group)
	.receiveBus_(channels[0][0])
	.delayTime_(0.1)
	.maxDelayTime_(0.5)
	.decayTime_(0.9)
	.sendBus_(0)
);

// #### INTRO NOISE
buffer = buffersC.at("plastic_ruffling_bag_1");
start = 0.0;
dur = 40;

score.add(
	playBufSynthDef.note(starttime: start, duration: dur, addAction: \head, target: group)
	.buffer_(buffer)
	.loop_(1)
	.rate_(1)
	.sendBus1_(channels[1][0])
	.sendBus2_(garbage)
	.gain_(10)
	.ampEnv_(CtkControl.env(Env([0, 1, 0.1, 0], [0.5, 0.01, 0.49], 6), starttime: start, timeScale: dur))
);

score.add(
	sosLPFilter.note(starttime: start, duration: dur, addAction: \tail, target: group)
	.freq_(CtkControl.env(Env([500, sampleRate / 2 - 10, 500], [0.5, 0.5], \lin), starttime: start, timeScale: dur))
	.receiveBus_(channels[1][0])
	.sendBus_(0)
);


// #### SHATTER AFTER VOCALS

buffer = buffersD.at("metal_hit_cage_1");
start = 19.8;
dur = buffer.duration * 5;

score.add(
	sgsSynth.note(starttime: start, duration: dur, addAction: \head, target: group)
	.dur_(dur)
	.gain_(12.0)
	.bufL_(getLeft.value(buffer))
	.bufR_(getRight.value(buffer))
	.wavRatio_(1.0)
	.sendBus1_(channels[0][0])
	.sendBus2_(0)
);

// ### extended shatter

buffer = buffersD.at("metal_hit_cage_1");
bufL = getLeft.value(buffer);
bufR = getRight.value(buffer);
start = 19.8 + 2.0;
dur = buffer.duration * 5 * 0.3 * 2;
loops = 30;

loops.do({arg i;
	score.add(
		sgsSynth.note(starttime: start + (i * dur - (dur * 0.1)), duration: dur, addAction: \head, target: group)
		.dur_(dur)
		.gain_(12.0)
		.startPos_(0.5)
		.endPos_(0.8)
		.bufL_(bufL)
		.bufR_(bufR)
		.wavRatio_(1.0)
		.sendBus1_(channels[0][0])
		.sendBus2_(garbage)
)});

// add multiple band passes
dur = dur * loops;

10.do({arg i;

	var freq =  400 * (i + 1);
	score.add(
		sosBPFilter.note(starttime: start, duration: dur, addAction: \tail, target: group)

		.freq_(CtkControl.env(Env.sinLFO(3, freq, freq + 300, 1, \sin), starttime: start, timeScale: dur))
		.bWidth_(800)
		.receiveBus_(channels[0][0])
		.sendBus_(channels[0][1])
	);
});

score.add(
	crossFade.note(starttime: start, duration: dur, addAction: \tail, target: group)
	.fade_(CtkControl.env(Env([-1, 1], [0.2], \lin), starttime: start, timeScale: dur))
	.receiveBusA_(channels[0][0])
	.receiveBusB_(channels[0][1])
	.sendBusA_(channels[0][2])
	.sendBusB_(channels[0][2])
);

score.add(
	combSpreader.note(starttime: start, duration: dur, addAction: \tail, target: group)
	.freq_(CtkControl.env(Env.sinLFO(6, 330, 500, 1, \sin), starttime: start, timeScale: dur))
	.b_(CtkControl.env(Env.sinLFO(6, 0.3, 0.7, 1, \sin), starttime: start, timeScale: dur))
	.widthAngle_(CtkControl.env(Env.sinLFO(7, 10, 30, 1, \sin), starttime: start, timeScale: dur))
	.receiveBus_(channels[0][2])
	.sendBus_(channels[0][3])
);

score.add(
	sosLPFilter.note(starttime: start, duration: dur, addAction: \tail, target: group)
	.freq_(CtkControl.env(Env([12000, 400], [1], \lin), starttime: start, timeScale: dur))
	.receiveBus_(channels[0][3])
	.sendBus_(0)
);

// ### PATTER OVER SHATTER
buffer = buffersC.at("paper_ripping_pages_1");
bufL = getLeft.value(buffer);
bufR = getRight.value(buffer);
start = 20;
dur = buffer.duration;

score.add(
	agsSynth.note(starttime: start, duration: dur, addAction: \head, target: group)
	.bufL_(bufL)
	.bufR_(bufR)
	.dur_(dur)
	.gain_(5.0)
	.sendBus1_(channels[2][0])
	.sendBus2_(0)
);

buffer = buffersC.at("rubber_pattering_balloon_1");
bufL = getLeft.value(buffer);
bufR = getRight.value(buffer);
start = 20;
dur = buffer.duration;
loops = 5;

loops.do({arg i;
	var thisStart = start + (dur * i);
	score.add(
		agsSynth.note(starttime: thisStart, duration: dur, addAction: \head, target: group)
		.bufL_(bufL)
		.bufR_(bufR)
		.dur_(dur)
		.gain_(-5.0)
		.sendBus1_(channels[2][0])
		.sendBus2_(garbage)
	);
});

dur = dur * loops;

score.add(
	vibTrem.note(starttime: start, duration: dur, addAction: \tail, target: group)
	.receiveBus_(channels[2][0])
	.modIndexGain_(-20)
	.rate_(CtkControl.env(Env([20, 1], [1], \lin), starttime: start, timeScale: dur))
	.minRate_(CtkControl.env(Env([20, 1], [1], \lin), starttime: start, timeScale: dur))
	.ratio_(CtkControl.env(Env([1.05, 1.01], [1], \lin), starttime: start, timeScale: dur))
	.sendBus_(channels[2][1])
);

score.add(
	universalComb.note(starttime: start, duration: dur, addAction: \tail, target: group)
	.freq_(CtkControl.env(Env([1000, 100], [1], \exp), starttime: start, timeScale: dur))
	.minFreq_(20)
	.b_(CtkControl.env(Env([0.8, 0.3], [1], \exp), starttime: start, timeScale: dur))
	.receiveBus_(channels[2][1])
	.sendBus_(0)
);

// ## MOOOOOTIONNNNN
buffer = buffersD.at("vocal_talking_motion_1");
start = 50.2;
dur = 4.0;
score.add(
	sgsSynth.note(starttime: start, duration: dur, addAction: \tail, target: group)
	.dur_(dur)
	.wavRatio_(1.0)
	.refFreq_(100)
	.freq_(100)
	.gain_(12.0)
	.startPos_(0.91)
	.endPos_(1.0)
	.bufL_(getLeft.value(buffer))
	.bufR_(getRight.value(buffer))
	.sendBus1_(channels[3][0])
	.sendBus2_(0)
	.ampEnv_(CtkControl.env(basicAmpEnv, starttime: start, timeScale: dur))
);

// reverb
score.add(
	dReverb.note(starttime: start, duration: 7.0, addAction: \tail, target: group)
	.filter_(CtkBuffer.collection(lopass.value(1024, 500, sampleRate) + hipass.value(1024, 5000, sampleRate)).addTo(score))
	.lowCut_(sampleRate / 2)
	.preDelayTime_(0.1)
	.maxDelayTime_(0.5)
	.receiveBus_(channels[3][0])
	.sendBus_(0)
);


// ## HIT NUMBER TWO


buffer = buffersD.at("water_plucking_bottle_1");
start = 56.3;
dur = 0.4 * 4;

score.add(
	playBufSynthDef.note(starttime: start, duration: dur, addAction: \head, target: group)
	.buffer_(buffer)
	.gain_(12.0)
	.startPos_(1.97)
	.rate_(0.25)
	.sendBus1_(channels[1][0])
	.sendBus2_(garbage)
);

score.add(
	disShaper.note(starttime: start, duration: dur, addAction: \tail, target: group)
	.receiveBus_(channels[1][0])
	.sendBus_(channels[1][2])
);

score.add(
	feedbackEcho.note(starttime: start, duration: dur + 5.0, addAction: \tail, target: group)
	.receiveBus_(channels[1][2])
	.sendBus_(0)
	.delayTime_(0.05)
	.maxDelayTime_(0.5)
	.decayTime_(2.0)
);


// ## metal ball rolling here

buffer = buffersC.at("metal_rolling_ball_3");
bufL = getLeft.value(buffer);
bufR = getRight.value(buffer);
start = 56.3;
dur = buffer.duration;

score.add(
	sgsSynth.note(starttime: start, duration: dur, addAction: \tail, target: group)
	.dur_(dur)
	.gain_(5.0)
	.bufL_(bufL)
	.bufR_(bufR)
	.wavRatio_(CtkControl.env(Env.sinLFO(20, Env([1, 0.1], [1], \lin), Env([1, 2], [1], \lin), 1.0, \sin), starttime: start, timeScale: dur))
	.sendBus1_(channels[2][0])
	.sendBus2_(0)
);

// add echo for the first 1.6 seconds with second hit
score.add(
	feedbackEcho.note(starttime: start, duration: 0.4 * 4, addAction: \tail, target: group)
	.receiveBus_(channels[2][0])
	.sendBus_(0)
	.delayTime_(0.05)
	.maxDelayTime_(0.5)
	.decayTime_(2.0)
);

score.add(
	weaverSSB.note(starttime: start, duration: dur, addAction: \tail, target: group)
	.receiveBus_(channels[2][0])
	.sendBus_(0)
	.ssbFreq_(CtkControl.env(Env([1.0, 0.1], [1], \exp), starttime: start, timeScale: dur))
	.bw_(CtkControl.env(Env.sinLFO(5, 0.1, 0.9, 0.5, \exp), starttime: start, timeScale: dur))
);

score.add(
	rotarySpreader.note(starttime: start, duration: dur, addAction: \tail, target: group)
	.receiveBus_(channels[2][0])
	.sendBus_(0)
	.ratio_(CtkControl.env(Env.sinLFO(5, 0.1, 0.9, 0.5, \exp), starttime: start, timeScale: dur))
	.maxRatio_(CtkControl.env(Env.sinLFO(5, 0.1, 0.9, 0.5, \exp), starttime: start, timeScale: dur))
	.rate_(CtkControl.env(Env([20, 2], [1.0], \lin), starttime: start, timeScale: dur))
	.minRate_(CtkControl.env(Env([10, 1], [1.0], \lin), starttime: start, timeScale: dur))
	.changeRate_(15)
	.change_(0.1)
	.rotateAngle_(CtkControl.env(Env.sinLFO(3, -45, 45, 1.0, \sin), starttime: start, timeScale: dur))
);

// ## shimmering dish stuff

buffer = buffersC.at("metal_shimmering_dish_1");
start = 99;
dur = buffer.duration;

score.add(
	playBufSynthDef.note(starttime: start, duration: dur, addAction: \tail, target: group)
	.dur_(dur)
	.buffer_(buffer)
	.sendBus1_(channels[2][0])
	.sendBus2_(garbage)
);

score.add(
	rotarySpreader.note(starttime: start, duration: dur, addAction: \tail, target: group)
	.receiveBus_(channels[2][0])
	.sendBus_(0)
	.rate_(CtkControl.env(Env([20, 3.0], [1], \exp), starttime: start, timeScale: dur))
	.minRate_(CtkControl.env(Env([20, 3.0], [1], \exp), starttime: start, timeScale: dur))
	.ratio_(0.8)
	.maxRatio_(0.8)
);

score.add(
	feedbackEcho.note(starttime: start + 3.0, duration: 10.0, addAction: \tail, target: group)
	.receiveBus_(channels[2][0])
	.sendBus_(0)
	.delayTime_(0.1)
	.maxDelayTime_(0.5)
	.decayTime_(3.0)
);


// ### RESONATING BOWL (under neath the metal ball)

buffer = buffersC.at("metal_resonating_bowl_1");
bufL = getLeft.value(buffer);
bufR = getRight.value(buffer);
start = 57;
dur = buffer.duration;

score.add(
	sgsSynth.note(starttime: start, duration: dur, addAction: \tail, target: group)
	.dur_(dur)
	.bufL_(bufL)
	.bufR_(bufR)
	.sendBus1_(channels[1][0])
	.sendBus2_(channels[1][1])
	.wavRatio_(CtkControl.env(Env.sinLFO(Env([1, 20, 1], [0.5, 0.5], \lin), 0.5, 1, 1.0, \sin), starttime: start, timeScale: dur))
);

score.add(
	firFilter.note(starttime: start, duration: dur, addAction: \tail, target: group)
	.receiveBus_(channels[1][0])
	.sendBus_(channels[1][2])
	.buffer_(CtkBuffer.collection(loshelf.value(1024, 1500, 3.0, sampleRate)).addTo(score))
);

// other intermediate thing here
score.add(
	sosBPFilter.note(starttime: start, duration: dur, addAction: \tail, target: group)
	.receiveBus_(channels[1][1])
	.sendBus_(channels[1][2])
	.freq_(CtkControl.env(Env.sinLFO(Env([3, 30, 3], [0.5, 0.5], \exp), 1500, 15000, 1, \exp), starttime: start, timeScale: dur))
	.bWidth_(400)
);

score.add(
	widthImage.note(starttime: start, duration: dur, addAction: \tail, target: group)
	.receiveBus_(channels[1][2])
	.sendBus_(0)
	.angle_(CtkControl.env(Env.sinLFO(Env([3, 30, 3], [0.5, 0.5], \exp), 10, 45, 1, \exp), starttime: start, timeScale: dur))
);


/// ######### PART 2

buffer = buffersD.at("sound_off_1");
bufL = getLeft.value(buffer);
bufR = getRight.value(buffer);
start = 107;
dur = buffer.duration;

score.add(
	sgsSynth.note(starttime: start, duration: dur, addAction: \tail, target: group)
	.dur_(dur)
	.gain_(8.0)
	.bufL_(bufL)
	.bufR_(bufR)
	.wavRatio_(CtkControl.env(Env([0.1, 1.0], [0.2], \exp), starttime: start, timeScale: dur))
	.sendBus1_(channels[3][0])
	.sendBus2_(0)
);

score.add(
	dReverb.note(starttime: start, duration: dur + 2.0, addAction: \tail, target: group)
	.receiveBus_(channels[3][0])
	.sendBus_(0)
	.lowCut_(CtkControl.env(Env([200, 600], [1.0], \exp), starttime: start, timeScale: dur))
	.filter_(allPassBuffer)
	.preDelayTime_(CtkControl.env(Env([0, 0.2], [1.0], 5), starttime: start, timeScale: dur))
	.maxDelayTime_(CtkControl.env(Env([0.1, 0.5], [1.0], \exp), starttime: start, timeScale: dur))
);


// MAKE SOME BRASSAGE


// ### subtle beeping that builds and terifies souls
buffer = buffersD.at("beep_1");
start = 120;
loops = 60;
dur = buffer.duration * loops;

// regular speed

score.add(
	playBufSynthDef.note(starttime: start, duration: dur, addAction: \tail, target: group)
	.buffer_(buffer)
	.dur_(dur)
	.gain_(8.0)
	.loop_(1)
	.rate_(1)
	.sendBus1_(channels[0][0])
	.sendBus2_(garbage)
);

score.add(
	sosHPFilter.note(starttime: start, duration: dur, addAction: \tail, target: group)
	.receiveBus_(channels[0][0])
	.sendBus_(0)
	.freq_(400)
);

score.add(
	dReverb.note(starttime: start, duration: dur + 2.0, addAction: \tail, target: group)
	.receiveBus_(channels[0][1])
	.sendBus_(0)
	.lowCut_(CtkControl.env(Env.noise1LFO(Env([1, 20], [1], \exp), 100, 1000, 1, \sin), starttime: start, timeScale: dur))
	.preDelayTime_(CtkControl.env(Env.noise1LFO(Env([1, 20], [1], \exp), 0.1, 0.3, 1, \sin), starttime: start, timeScale: dur))
	.filter_(allPassBuffer)
);



// two octaves below

score.add(
	playBufSynthDef.note(starttime: start + 10, duration: dur - 10, addAction: \tail, target: group)
	.buffer_(buffer)
	.dur_(dur)
	.gain_(8.0)
	.loop_(1)
	.rate_(0.25)
	.sendBus1_(channels[1][0])
	.sendBus2_(0)
);

score.add(
	dReverb.note(starttime: start, duration: dur + 2.0, addAction: \tail, target: group)
	.receiveBus_(channels[1][0])
	.sendBus_(0)
	.lowCut_(CtkControl.env(Env.noise1LFO(Env([1, 20], [1], \exp), 100, 1000, 1, \sin), starttime: start, timeScale: dur))
	.filter_(allPassBuffer)
);


// one octave higher

score.add(
	playBufSynthDef.note(starttime: start + 30, duration: dur - 30, addAction: \tail, target: group)
	.buffer_(buffer)
	.dur_(dur)
	.gain_(3.0)
	.loop_(1)
	.rate_(2.0)
	.sendBus1_(channels[2][0])
	.sendBus2_(0)
);

// abs beeping

score.add(
	absShaper.note(starttime: start + 40, duration: dur - 40, addAction: \tail, target: group)
	.receiveBus_(channels[1][0])
	.sendBus_(0)
);

// echo out of all

score.add(
	feedbackEcho.note(starttime: start + 55, duration: 10.0, addAction: \tail, target: group)
	.receiveBus_(channels[1][0])
	.sendBus_(0)
	.decayTime_(2.0)
	.delayTime_(0.1)
	.maxDelayTime_(0.15)
);

score.add(
	feedbackEcho.note(starttime: start + 58, duration: 10.0, addAction: \tail, target: group)
	.receiveBus_(channels[0][0])
	.sendBus_(0)
	.decayTime_(5.0)
	.delayTime_(0.2)
	.maxDelayTime_(0.3)
);


// ## vocals are back
buffer = buffersD.at("vocal_talking_motion_1");
start = 120;
dur = 60;

score.add(
	sgsSynth.note(starttime: start, duration: dur, addAction: \head, target: group)
	.dur_(dur)
	.wavRatio_(CtkControl.env(Env.sinLFO(10, 0.5, 1, 1, \sin), starttime: start, timeScale: dur))
	.refFreq_(100)
	.freq_(CtkControl.env(Env.sinLFO(10, 50, 100, 1, \sin), starttime: start, timeScale: dur))
	.gain_(12.0)
	.bufL_(getLeft.value(buffer))
	.bufR_(getRight.value(buffer))
	.sendBus1_(channels[4][0])
	.sendBus2_(garbage)
);

score.add(
	sosLPFilter.note(starttime: start, duration: dur, addAction: \tail, target: group)
	.receiveBus_(channels[4][0])
	.sendBus_(0)
	.freq_(CtkControl.env(Env.sinLFO(20, 100, 300, 1, \sin), starttime: start, timeScale: dur))
);

// high pass filtered

score.add(
	sgsSynth.note(starttime: start + 13, duration: dur - 13, addAction: \head, target: group)
	.dur_(dur)
	.wavRatio_(CtkControl.env(Env.sinLFO(10, 1, 2.0, 1, \sin), starttime: start, timeScale: dur))
	.refFreq_(100)
	.freq_(CtkControl.env(Env.sinLFO(10, 100, 200, 1, \sin), starttime: start, timeScale: dur))
	.gain_(12.0)
	.bufL_(getLeft.value(buffer))
	.bufR_(getRight.value(buffer))
	.sendBus1_(channels[3][0])
	.sendBus2_(garbage)
);

score.add(
	sosHPFilter.note(starttime: start + 13, duration: dur - 13, addAction: \tail, target: group)
	.receiveBus_(channels[3][0])
	.sendBus_(0)
	.freq_(CtkControl.env(Env.sinLFO(Env([20, 1], [1], \sin), 6000, 7000, 1, \sin), starttime: start, timeScale: dur))
);

score.add(
	sosHPFilter.note(starttime: start + 13, duration: dur - 13, addAction: \tail, target: group)
	.receiveBus_(channels[3][0])
	.sendBus_(0)
	.freq_(CtkControl.env(Env.sinLFO(Env([1, 20], [1], \sin), 8000, 10000, 1, \sin), starttime: start, timeScale: dur))
);

// distorted

score.add(
	sgsSynth.note(starttime: start + 30, duration: buffer.duration + 5.0, addAction: \tail, target: group)
	.dur_(dur)
	.wavRatio_(CtkControl.env(Env.sinLFO(10, 0.9, 1, 1, \sin), starttime: start, timeScale: dur))
	.refFreq_(100)
	.freq_(CtkControl.env(Env.sinLFO(10, 50, 100, 1, \sin), starttime: start, timeScale: dur))
	.gain_(14.0)
	.bufL_(getLeft.value(buffer))
	.bufR_(getRight.value(buffer))
	.sendBus1_(channels[1][8])
	.sendBus2_(garbage)
);

score.add(
	cosShaper.note(starttime: start + 30, duration: buffer.duration, addAction: \tail, target: group)
	.receiveBus_(channels[1][8])
	.sendBus_(0)
);

// some nice combing
buffer = buffersC.at("paper_crumpling_wrapper_1");
start = 120;
dur = buffer.duration * 2;

score.add(
	playBufSynthDef.note(starttime: start, duration: dur, addAction: \tail, target: group)
	.gain_(-3.0)
	.dur_(dur)
	.loop_(1)
	.buffer_(buffer)
	.rate_(1)
	.sendBus1_(garbage)
	.sendBus2_(channels[1][5])
	.ampEnv_(CtkControl.env(basicAmpEnv, starttime: start, timeScale: dur))
);

score.add(
	combSpreader.note(starttime: start, duration: dur, addAction: \tail, target: group)
	.evenLeft_(0)
	.freq_(CtkControl.env(Env.sinLFO(2, 400, 733, 1, \sin), starttime: start, timeScale: dur))
	.minFreq_(440)
	.b_(CtkControl.env(Env.sinLFO(1, 0.1, 0.7, 1, \sin), starttime: start, timeScale: dur))
	.receiveBus_(channels[1][5])
	.sendBus_(0)
	.widthAngle_(45)
);

// finish with "YOU CAN'T TURN THAT SOUND OFF"

buffer = buffersD.at("sound_off_1");
start = 181;
dur = (1.5 / 0.95);
score.add(
	playBufSynthDef.note(starttime: start, duration: dur, addAction: \tail, target: group)
	.dur_(dur)
	.startPos_(7.84)
	.buffer_(buffer)
	.rate_(0.95)
	.loop_(0)
	.sendBus1_(channels[2][6])
	.sendBus2_(garbage)
);

score.add(
	firFilter.note(starttime: start, duration: dur, addAction: \tail, target: group)
	.receiveBus_(channels[2][6])
	.sendBus_(channels[2][7])
	.buffer_(CtkBuffer.collection(lopass.value(1024, 6000, sampleRate)).addTo(score))
);

score.add(
	widthImage.note(starttime: start, duration: dur, addAction: \tail, target: group)
	.receiveBus_(channels[2][7])
	.sendBus_(0)
	.angle_(CtkControl.env(Env([45, 0], [1], \exp), starttime: start, timeScale: dur))
);

score.add(
	dReverb.note(starttime: start, duration: dur + 1.0, addAction: \tail, target: group)
	.receiveBus_(channels[2][7])
	.sendBus_(0)
	.filter_(allPassBuffer)
);

///////////////// RENDER THE SCORE //////////////////
// write score to sound file with the -write message
// NOTE: we're using argument names to specify the args. For 'duration', we're letting Ctk
//       do the work for us!
score.write(
	path: outputPath.standardizePath,
	sampleRate: sampleRate,
	headerFormat: headerFormat,
	sampleFormat: sampleFormat,
	options: ServerOptions.new.numOutputBusChannels_(numOutputChannels)
);

///////////////// FREE BUSSES //////////////////
channels.asArray.do({arg item, i;
	item.free
});
garbage.free;
)
SFPlayer("Gaimari_final.wav".resolveRelative).gui;