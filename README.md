# spindle-ripple-coupling-
determination of TFR(time-frequency representations) and PAC ( Phase amplitude coupling ) and SI (synchronization index  )

I used Thalamus and hippocampus's data 
Firstly, I separated the seconds of each event from rem, nrem, wake.Then, the data is filtered by 7-16 HZ and fs=1250 frequency.I showed the parts of NREM,REM,WAKE with numbers 3, 2 and 1,respectively.
I used the range of 0.5 to 300 Hz with 3 cycles of the low frequency cut off and finite-impulse-response (FIR)  to display the signal.I separated the data and times of NREM for both channels (Thalamus and hippocampus) and then filtered them between 7-16 HZ And using function FINDSPINLDE detected spindles.
subsequently, I separated NREM data which their stop times were less than start times of REM data and I calculate spindles of  25 S before REM's start .finally I calculated their average.

# The method to calculate phase amplitude coupling (PAC):
 First time–frequency representations (TFRs) were calculated for every spindle event by mtmconvol function of the FieldTrip toolbox with frequency steps of 0.25 Hz in the frequency range of 5-300 Hz. Sliding (10 ms steps) Hanning  
tapered  window with a variable length including five cycles were used to ensure reliable power estimates. TFRs of all epochs around an event were then normalized as percentage change from  the pre-event baseline (−2.5 to −1.5 s). 
To quantify the modulation of hippocampal ripple amplitude with the phase of thalamic spindle, first TFR bins around each spindle were averaged across the ripple frequency ranges to obtain the power of ripple time series (-3 to 3 s around the spindle peak). To ensure proper phase estimation, both Thalamus and hippocampus were filtered in the spindle frequency range.. Next,we extracted the phase values of these time series using the Hilbert transform and defined a synchronization index (SI) for each spindle event as the vector mean of phase difference between the two time series.

# The method to calculate SI :
![Capture](https://github.com/mohamad9014/spindle-ripple-coupling-/assets/121359931/e4509862-f7b6-49ff-bc50-5506d0522fec)

where m is the number of time points, θ𝑠𝑝𝑖𝑛𝑑𝑙𝑒(j)  represents the phase value of the spindle time series at  time point j and θ𝑟𝑖𝑝𝑝𝑙𝑒(j) shows the phase value of the fluctuations in the ripple power time series at time point j.
The used interval for estimating the preferred phase was −0.25 s to +0.25 s around the spindle peak. The resulting SI is a complex number with the angle representing  the phase shift between the two oscillatory events (i.e. ripple power and spindle) while the absolute value of SI represents the strength of the coupling between the two oscillatory events. It gives me p-value and z-value .If its p-value is less than 0.05 my data is oriented .

# TFR figure :
![Capture-2](https://github.com/mohamad9014/spindle-ripple-coupling-/assets/121359931/bd6bb2dc-88db-4af7-976c-b56b06e35d27)

As can be seen in the figure, 0 degrees is the peak and 180 degrees is the valley And we have a very strong power (ripplr coupling ) before the peak (zero degree)  , which corresponds to the figure below:

![Capture-3](https://github.com/mohamad9014/spindle-ripple-coupling-/assets/121359931/c63e2d1b-8603-436e-a4a0-29c926135ab8)


