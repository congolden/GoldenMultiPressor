#pragma once

#ifndef __fourWaySplitter__
#define __fourWaySplitter__

#include "fxobjects.h"
#include "math.h"

/**
\struct fourWaySplitterParameters
\ingroup FX-Objects
\brief
Custom parameter structure for the fourWaySplitter object.

\author <Your Name> <http://www.yourwebsite.com>
\remark <Put any remarks or notes here>
\version Revision : 1.0
\date Date : 2019 / 01 / 31
*/
struct fourWaySplitterParameters
{
	fourWaySplitterParameters() {}

	/** all FXObjects parameter objects require overloaded= operator so remember to add new entries if you add new variables. */
	fourWaySplitterParameters& operator=(const fourWaySplitterParameters& params)	// need this override for collections to work
	{
		// --- it is possible to try to make the object equal to itself
		//     e.g. thisObject = thisObject; so this code catches that
		//     trivial case and just returns this object
		if (this == &params)
			return *this;

		// --- copy from params (argument) INTO our variables
		//myVariable = params.myVariable;
		L_freq = params.L_freq;
		M_freq = params.M_freq;
		H_freq = params.H_freq;

		inputVolume_dB = params.inputVolume_dB;
		outputVolume_dB = params.outputVolume_dB;
		midsideEnable = params.midsideEnable;
		sideSelect = params.sideSelect;

		LF_enable = params.LF_enable;
		LF_threshold_dB= params.LF_threshold_dB;
		LF_ratio = params.LF_ratio;
		LF_attackTime_ms = params.LF_attackTime_ms;
		LF_releaseTime_ms = params.LF_releaseTime_ms;
		LF_makeupGain_dB = params.LF_makeupGain_dB;
		LF_soloEnable = params.LF_soloEnable;
		LF_grMeter = params.LF_grMeter;
		LFL_volMeter = params.LFL_volMeter;
		LFR_volMeter = params.LFR_volMeter;

		LM_enable = params.LM_enable;
		LM_threshold_dB = params.LM_threshold_dB;
		LM_ratio = params.LM_ratio;
		LM_attackTime_ms = params.LM_attackTime_ms;
		LM_releaseTime_ms = params.LM_releaseTime_ms;
		LM_makeupGain_dB = params.LM_makeupGain_dB;
		LM_soloEnable = params.LM_soloEnable;
		LM_grMeter = params.LM_grMeter;
		LML_volMeter = params.LML_volMeter;
		LMR_volMeter = params.LMR_volMeter;


		HM_enable = params.HM_enable;
		HM_threshold_dB = params.HM_threshold_dB;
		HM_ratio = params.HM_ratio;
		HM_attackTime_ms = params.HM_attackTime_ms;
		HM_releaseTime_ms = params.HM_releaseTime_ms;
		HM_makeupGain_dB = params.HM_makeupGain_dB;
		HM_soloEnable = params.HM_soloEnable;
		HM_grMeter = params.HM_grMeter;
		HML_volMeter = params.HML_volMeter;
		HMR_volMeter = params.HMR_volMeter;

		HF_enable = params.HF_enable;
		HF_threshold_dB = params.HF_threshold_dB;
		HF_ratio = params.HF_ratio;
		HF_attackTime_ms = params.HF_attackTime_ms;
		HF_releaseTime_ms = params.HF_releaseTime_ms;
		HF_makeupGain_dB = params.HF_makeupGain_dB;
		HF_soloEnable = params.HF_soloEnable;
		HF_grMeter = params.HF_grMeter;
		HFL_volMeter = params.HFL_volMeter;
		HFR_volMeter = params.HFR_volMeter;

		// --- MUST be last
		return *this;
	}
	double L_freq = 200;
	double M_freq = 1000;
	double H_freq = 8000;

	double inputVolume_dB = 0.0;
	double outputVolume_dB = 0.0;

	double LF_enable = true;
	double LF_threshold_dB = 0.0;
	double LF_ratio = 0.0;
	double LF_attackTime_ms = 0.0;
	double LF_releaseTime_ms = 0.0;
	double LF_makeupGain_dB = 0.0;
	double LF_soloEnable = false;
	float LF_grMeter = 0.f;
	float LFL_volMeter = 0.f;
	float LFR_volMeter = 0.f;

	double LM_enable = true;
	double LM_threshold_dB = 0.0;
	double LM_ratio = 0.0;
	double LM_attackTime_ms = 0.0;
	double LM_releaseTime_ms = 0.0;
	double LM_makeupGain_dB = 0.0;
	double LM_soloEnable = false;
	float LM_grMeter = 0.f;
	float LML_volMeter = 0.f;
	float LMR_volMeter = 0.f;

	double HM_enable = true;
	double HM_threshold_dB = 0.0;
	double HM_ratio = 0.0;
	double HM_attackTime_ms = 0.0;
	double HM_releaseTime_ms = 0.0;
	double HM_makeupGain_dB = 0.0;
	double HM_soloEnable = false;
	float HM_grMeter = 0.f;
	float HML_volMeter = 0.f;
	float HMR_volMeter = 0.f;

	double HF_enable = true;
	double HF_threshold_dB = 0.0;
	double HF_ratio = 0.0;
	double HF_attackTime_ms = 0.0;
	double HF_releaseTime_ms = 0.0;
	double HF_makeupGain_dB = 0.0;
	double HF_soloEnable = false;
	float HF_grMeter = 0.f;
	float HFL_volMeter = 0.f;
	float HFR_volMeter = 0.f;

	bool midsideEnable = false;
	bool sideSelect = false;
	///////

	// --- individual parameters
	//data_type myVariable = 0.0;	///< init
};


/**
\class fourWaySplitter
\ingroup FX-Objects
\brief
The fourWaySplitter object implements ....

Audio I/O:
- Processes mono input to mono output.
- *** Optionally, process frame *** Modify this according to your object functionality

Control I/F:
- Use fourWaySplitterParameters structure to get/set object params.

\author <Your Name> <http://www.yourwebsite.com>
\remark <Put any remarks or notes here>
\version Revision : 1.0
\date Date : 2019 / 01 / 31
*/
class fourWaySplitter : public IAudioSignalProcessor
{
public:
	/* C-TOR */
	fourWaySplitter()
	{
		/*	AudioFilterParameters params = lpFilter.getParameters();
		params.algorithm = filterAlgorithm::kLWRLPF2;
		lpFilter.setParameters(params);

		params = lmhpFilter.getParameters();
		params.algorithm = filterAlgorithm::kLWRHPF2;
		lmhpFilter.setParameters(params);

		params = lmlpFilter.getParameters();
		params.algorithm = filterAlgorithm::kLWRLPF2;
		lmlpFilter.setParameters(params);

		params = hmhpFilter.getParameters();
		params.algorithm = filterAlgorithm::kLWRHPF2;
		hmhpFilter.setParameters(params);

		params = hmlpFilter.getParameters();
		params.algorithm = filterAlgorithm::kLWRLPF2;
		hmlpFilter.setParameters(params);

		params = hpFilter.getParameters();
		params.algorithm = filterAlgorithm::kLWRHPF2;
		hpFilter.setParameters(params);
		*/
	}	
	~fourWaySplitter(void) {}	/* D-TOR */

public:
	/** reset members to initialized state */
	virtual bool reset(double _sampleRate)
	{
		// --- store the sample rate
		sampleRate = _sampleRate;

		splitterFilters->reset(_sampleRate);

		LFL_compressor.reset(_sampleRate);
		LFR_compressor.reset(_sampleRate);
		LML_compressor.reset(_sampleRate);
		LMR_compressor.reset(_sampleRate);
		HML_compressor.reset(_sampleRate);
		HMR_compressor.reset(_sampleRate);
		HFL_compressor.reset(_sampleRate);
		HFR_compressor.reset(_sampleRate);
		// --- do any other per-audio-run inits here

		return true;
	}


	/** process MONO input */
	/**
	\param xn input
	\return the processed sample
	*/
	virtual double processAudioSample(double xn)
	{
		// --- the output variable
		double yn = 0.0;

		// --- do your DSP magic here to create yn


		// --- done
		return yn;
	}

	/** query to see if this object can process frames */
	virtual bool canProcessAudioFrame() { return true; } // <-- change this!

	/** process audio frame: implement this function if you answer "true" to above query */
	virtual bool processAudioFrame(const float* inputFrame,	/* ptr to one frame of data: pInputFrame[0] = left, pInputFrame[1] = right, etc...*/
					     float* outputFrame,
					     uint32_t inputChannels,
					     uint32_t outputChannels)


	{
		double xnL = inputFrame[0];
		double xnR = inputFrame[1];

		double mid;
		double sides;

		double ynL = xnL * inputVolume_cooked;
		double ynR = xnR * inputVolume_cooked;

		 // mid side encode

		if (parameters.midsideEnable)
		{
			mid = .5*(ynL+ynR);
			sides = .5*(ynL-ynR);

			ynL = mid;
			ynR = sides;
		}
		

		FilterBankOutput HPF_Left = splitterFilters[10].processFilterBank(ynL);
		FilterBankOutput HPF_Right = splitterFilters[11].processFilterBank(ynR);
		
		
		FilterBankOutput LPF_Left = splitterFilters[0].processFilterBank(ynL);
		FilterBankOutput LPF_Right = splitterFilters[1].processFilterBank(ynR);


		FilterBankOutput LMF_Left = splitterFilters[4].processFilterBank(LPF_Left.HFOut);
		FilterBankOutput LMF_Right = splitterFilters[5].processFilterBank(LPF_Right.HFOut);

		FilterBankOutput HMF_Left = splitterFilters[6].processFilterBank(HPF_Left.LFOut);
		FilterBankOutput HMF_Right = splitterFilters[7].processFilterBank(HPF_Right.LFOut);

		double LF_L = LPF_Left.LFOut;
		double LF_R = LPF_Right.LFOut;
		double LM_L = LMF_Left.LFOut;
		double LM_R = LMF_Right.LFOut;
		double HM_L = (-HMF_Left.HFOut);
		double HM_R = (-HMF_Right.HFOut);
		double HF_L = HPF_Left.HFOut;
		double HF_R = HPF_Right.HFOut;

		// dynamics
		if (parameters.LF_enable)
		{
			DynamicsProcessorParameters LFL_params = LFL_compressor.getParameters();
			DynamicsProcessorParameters LFR_params = LFR_compressor.getParameters();
			if (parameters.midsideEnable)
			{
				if (parameters.sideSelect) // just compress side
				{
					LF_R = LFR_compressor.processAudioSample(LF_R);
					parameters.LF_grMeter = (-LFR_params.gainReduction + 1) / 2;
				}
				else   // just compress mid
				{
					LF_L = LFL_compressor.processAudioSample(LF_L);
					parameters.LF_grMeter = (-LFL_params.gainReduction + 1) / 2;
				}
			}
			else //compres both
			{
				LF_L = LFL_compressor.processAudioSample(LF_L);
				LF_R = LFR_compressor.processAudioSample(LF_R);
				parameters.LF_grMeter = (-LFL_params.gainReduction + 1) / 2;
			}
			
		}
		else
		{
			LF_L = 0;
			LF_R = 0;
		}
		parameters.LFL_volMeter = LF_L*2;
		parameters.LFR_volMeter = LF_R*2;

		if (parameters.LM_enable)
		{
			DynamicsProcessorParameters LML_params = LML_compressor.getParameters();
			DynamicsProcessorParameters LMR_params = LMR_compressor.getParameters();
			if (parameters.midsideEnable)
			{
				if (parameters.sideSelect) // just compress side
				{
					LM_R = LMR_compressor.processAudioSample(LM_R);
					parameters.LM_grMeter = (-LMR_params.gainReduction + 1) / 2;
				}
				else   // just compress mid
				{
					LM_L = LML_compressor.processAudioSample(LM_L);
					parameters.LM_grMeter = (-LML_params.gainReduction + 1) / 2;
				}
			}
			else //compres both
			{
				LM_L = LML_compressor.processAudioSample(LM_L);
				LM_R = LMR_compressor.processAudioSample(LM_R);
				parameters.LM_grMeter = (-LML_params.gainReduction + 1) / 2;
			}
		}
		else
		{
			LM_L = 0;
			LM_R = 0;
		}
		parameters.LML_volMeter = LM_L*2;
		parameters.LMR_volMeter = LM_R*2;

		if (parameters.HM_enable)
		{
			DynamicsProcessorParameters HML_params = HML_compressor.getParameters();
			DynamicsProcessorParameters HMR_params = HMR_compressor.getParameters();
			if (parameters.midsideEnable)
			{
				if (parameters.sideSelect) // just compress side
				{
					HM_R = HMR_compressor.processAudioSample(HM_R);
					parameters.HM_grMeter = (-HMR_params.gainReduction + 1) / 2;
				}
				else   // just compress mid
				{
					HM_L = HML_compressor.processAudioSample(HM_L);
					parameters.HM_grMeter = (-HML_params.gainReduction + 1) / 2;
				}
			}
			else //compres both
			{
				HM_L = HML_compressor.processAudioSample(HM_L);
				HM_R = HMR_compressor.processAudioSample(HM_R);
				parameters.HM_grMeter = (-HML_params.gainReduction + 1) / 2;
			}
		}
		else
		{
			HM_L = 0;
			HM_R = 0;
		}
		parameters.HML_volMeter = HM_L*2;
		parameters.HMR_volMeter = HM_R*2;

		if (parameters.HF_enable)
		{
			DynamicsProcessorParameters HFL_params = HFL_compressor.getParameters();
			DynamicsProcessorParameters HFR_params = HFR_compressor.getParameters();
			if (parameters.midsideEnable)
			{
				if (parameters.sideSelect) // just compress side
				{
					HF_R = HFR_compressor.processAudioSample(HF_R);
					parameters.HF_grMeter = (-HFR_params.gainReduction + 1) / 2;
				}
				else   // just compress mid
				{
					HF_L = HFL_compressor.processAudioSample(HF_L);
					parameters.HF_grMeter = (-HFL_params.gainReduction + 1) / 2;
				}
			}
			else //compres both
			{
				HF_L = HFL_compressor.processAudioSample(HF_L);
				HF_R = HFR_compressor.processAudioSample(HF_R);
				parameters.HF_grMeter = (-HFL_params.gainReduction + 1) / 2;
			}
		}
		else
		{
			HF_L = 0;
			HF_R = 0;
		}
		parameters.HFL_volMeter = HF_L*2;
		parameters.HFR_volMeter = HF_R*2;
		//Volume

		double sumL = LF_L + LM_L + HM_L + HF_L;  // or mid
		double sumR = LF_R + LM_R + HM_R + HF_R;  // or sides
		//mid side decode
		
		if (parameters.midsideEnable)
		{
			sumL = sumL + sumR;
			sumR = sumL - sumR;
		}
		
		//output

		ynL = sumL *outputVolume_cooked;
		ynR = sumR *outputVolume_cooked;

	


		// --- mono/mono; use existing function and done
		if (inputChannels == 1 && outputChannels == 1)
		{
			outputFrame[0] = ynL;
			return true; // --- processed!
		}


		// --- mono/stereo: pan ledt input channel to left + right outputs
		if (inputChannels == 1 && outputChannels == 2)
		{
			outputFrame[0] = ynL;
			outputFrame[1] = ynL;
			// --- outbound variables

			return true; // --- processed
		}

		// --- stereo/stereo: pan is actually a balance control
		else if (inputChannels == 2 && outputChannels == 2)
		{
			outputFrame[0] = ynL;
			outputFrame[1] = ynR;
			// --- outboud variables
			return true; // --- processed!
		}
		return false; // NOT handled
	}


	/** get parameters: note use of custom structure for passing param data */
	/**
	\return fourWaySplitterParameters custom data structure
	*/
	fourWaySplitterParameters getParameters()
	{
		return parameters;
	}

	/** set parameters: note use of custom structure for passing param data */
	/**
	\param fourWaySplitterParameters custom data structure
	*/
	void setParameters(const fourWaySplitterParameters& params)
	{
		// --- copy them; note you may choose to ignore certain items
		//     and copy the variables one at a time, or you may test
		//     to see if cook-able variables have changed; if not, then
		//     do not re-cook them as it just wastes CPU
		parameters = params;

		// --- update member objects
		LRFilterBankParameters LPF_params;
		LRFilterBankParameters MPF_params;
		LRFilterBankParameters HPF_params;
		LPF_params.splitFrequency = parameters.L_freq;
		MPF_params.splitFrequency = parameters.M_freq;
		HPF_params.splitFrequency = parameters.H_freq;

		splitterFilters[0].setParameters(LPF_params);
		splitterFilters[1].setParameters(LPF_params);
		splitterFilters[2].setParameters(LPF_params);
		splitterFilters[3].setParameters(LPF_params);
		splitterFilters[4].setParameters(MPF_params);
		splitterFilters[5].setParameters(MPF_params);
		splitterFilters[6].setParameters(MPF_params);
		splitterFilters[7].setParameters(MPF_params);
		splitterFilters[8].setParameters(HPF_params);
		splitterFilters[9].setParameters(HPF_params);
		splitterFilters[10].setParameters(HPF_params);
		splitterFilters[11].setParameters(HPF_params);

		DynamicsProcessorParameters LF_params = LFL_compressor.getParameters();
		LF_params.threshold_dB = parameters.LF_threshold_dB;
		LF_params.ratio = parameters.LF_ratio;
		LF_params.attackTime_mSec = parameters.LF_attackTime_ms;
		LF_params.releaseTime_mSec = parameters.LF_releaseTime_ms;
		LF_params.outputGain_dB = parameters.LF_makeupGain_dB;
		LF_params.calculation == dynamicsProcessorType::kCompressor;
		
		LFL_compressor.setParameters(LF_params);
		LFR_compressor.setParameters(LF_params);

		DynamicsProcessorParameters LM_params = LML_compressor.getParameters();
		LM_params.threshold_dB = parameters.LM_threshold_dB;
		LM_params.ratio = parameters.LM_ratio;
		LM_params.attackTime_mSec = parameters.LM_attackTime_ms;
		LM_params.releaseTime_mSec = parameters.LM_releaseTime_ms;
		LM_params.outputGain_dB = parameters.LM_makeupGain_dB;
		LM_params.calculation == dynamicsProcessorType::kCompressor;
		LML_compressor.setParameters(LM_params);
		LMR_compressor.setParameters(LM_params);

		DynamicsProcessorParameters HM_params = HML_compressor.getParameters();
		HM_params.threshold_dB = parameters.HM_threshold_dB;
		HM_params.ratio = parameters.HM_ratio;
		HM_params.attackTime_mSec = parameters.HM_attackTime_ms;
		HM_params.releaseTime_mSec = parameters.HM_releaseTime_ms;
		HM_params.outputGain_dB = parameters.HM_makeupGain_dB;
		HM_params.calculation == dynamicsProcessorType::kCompressor;
		HML_compressor.setParameters(HM_params);
		HMR_compressor.setParameters(HM_params);

		DynamicsProcessorParameters HF_params = HFL_compressor.getParameters();
		HF_params.threshold_dB = parameters.HF_threshold_dB;
		HF_params.ratio = parameters.HF_ratio;
		HF_params.attackTime_mSec = parameters.HF_attackTime_ms;
		HF_params.releaseTime_mSec = parameters.HF_releaseTime_ms;
		HF_params.outputGain_dB = parameters.HF_makeupGain_dB;
		HF_params.calculation == dynamicsProcessorType::kCompressor;
		HFL_compressor.setParameters(HF_params);
		HFR_compressor.setParameters(HF_params);

		if (parameters.inputVolume_dB == -60.0)
			inputVolume_cooked = 0.0;
		else
			inputVolume_cooked = pow(10.0, parameters.inputVolume_dB / 20.0);

		if (parameters.outputVolume_dB == -60.0)
			outputVolume_cooked = 0.0;
		else
			outputVolume_cooked = pow(10.0, parameters.outputVolume_dB / 20.0);


		

		// --- cook parameters here
	}

private:
	fourWaySplitterParameters parameters; ///< object parameters
	/* AudioFilter lpFilter; ///< low-band filter
	AudioFilter lmhpFilter;
	AudioFilter lmlpFilter;
	AudioFilter hmhpFilter;
	AudioFilter hmlpFilter;
	AudioFilter hpFilter; ///< high-band filter
	*/
	LRFilterBank splitterFilters[12]; // 2 per band 
	DynamicsProcessor LFL_compressor;
	DynamicsProcessor LFR_compressor;
	DynamicsProcessor LML_compressor;
	DynamicsProcessor LMR_compressor;
	DynamicsProcessor HML_compressor;
	DynamicsProcessor HMR_compressor;
	DynamicsProcessor HFL_compressor;
	DynamicsProcessor HFR_compressor;
	// --- local variables used by this object
	double sampleRate = 0.0;	///< sample rate
	double inputVolume_cooked = 1.0;
	double outputVolume_cooked = 1.0;
	double hmf_vol_cooked = 1.0;
	double hf_vol_cooked = 1.0;
};

#endif