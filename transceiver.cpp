#include <iostream>
#include <vector>
#include <thread>
#include <atomic>
#include <complex>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "portaudio.h"
#include <cstring>
#include <cctype>
#include <chrono>
#include <future>

//ALWAYS ENSURE THAT THE SAMPLERATE AND SAMPLE_RATE ARE IDENTICAL!!
using namespace std;
typedef unsigned char u8;
typedef unsigned short u16;
typedef unsigned long u32;
typedef unsigned long long u64;
// Defs from FFT code
int FREQS[8]= {697,770,852,941,1209,1336,1477,1633};
// Frequencies for DTMF tones
float frqU[10] = {1336, 1209, 1336, 1477, 1209, 1336, 1477, 1209, 1336, 1477};
float frqD[10] = {941,  697,  697,  697,  770,  770,  770,  852,  852,  852};
// Lower two together, upper two together
float V21FRQ[4]={980,1180,1650,1850};
int issender;
int echoSuppressed=0;
// Global variables
vector<u8> inputNumbers;
vector<double> echosuppressionOut = {};

// Struct to pass frequency pairs
struct FrequencyPair {
    float frequency1;
    float frequency2;
};
FrequencyPair frequencies = {0.0f, 0.0f};

#define IGNORE_ERRORS // enable this to skip termination for NACK

#define SAMPLERATE 44100
#define THETA 6.28318530717958647692528676655900  // Constant: 2*pi
#define PI 3.14159265358979323846264338327950

#define SAMPLE_RATE  44100
#define FRAMES_PER_BUFFER (512)
#define NUM_CHANNELS    (2)

#define DITHER_FLAG     (0)
#define WRITE_TO_FILE   (1)
#define SAMPLE_SIZE sizeof(float)
#define DURATION 0.04f  //DTMF symbol duration

/* Select sample format. */
#if 1
#define PA_SAMPLE_TYPE  paFloat32
typedef float SAMPLE;
#define SAMPLE_SILENCE  (0.0f)
#define PRINTF_S_FORMAT "%.8f"
#elif 1
#define PA_SAMPLE_TYPE  paInt16
typedef short SAMPLE;
#define SAMPLE_SILENCE  (0)
#define PRINTF_S_FORMAT "%d"
#elif 0
#define PA_SAMPLE_TYPE  paInt8
typedef char SAMPLE;
#define SAMPLE_SILENCE  (0)
#define PRINTF_S_FORMAT "%d"
#else
#define PA_SAMPLE_TYPE  paUInt8
typedef unsigned char SAMPLE;
#define SAMPLE_SILENCE  (128)
#define PRINTF_S_FORMAT "%d"
#endif
/*-------------------------------START DETCTION SECTION HERE-----------------------------------*/

typedef struct{
    int          frameIndex;  /* Index into sample array. */
    int          maxFrameIndex;
    SAMPLE      *recordedSamples;
}
paTestData;

//Compute the CT-FFT of the input.
vector<complex<double>> f(vector<complex<double>> y){
    int n=y.size();
    if(n==1){return y;}
    vector<complex<double>> gamma;
    if(n&(n-1)){printf("Error -> choose N=2^j\n");}
    else{
        n=n/2;
        vector<complex<double>> a;
        vector<complex<double>> b;
        vector<complex<double>> c;
        for(int i=0;i<n;i++){
            a.push_back(y[2*i]);
            b.push_back(y[1+2*i]);
        }
        a=f(a);
        b=f(b);
        complex<double> t(0,-PI/n);
        int cx=0;
        for(complex<double> k(0,0);k.real()<n;k.real(k.real()+1)){
            c.push_back( exp(t*k)*b[cx]   );//b[k] might by complex -> don't integrate into polar.
            cx++;
        }

        for(int i=0;i<n;i++){
            gamma.push_back(a[i]+c[i]);
        }
        for(int i=0;i<n;i++){
            gamma.push_back(a[i]-c[i]);
        }
    }
    return gamma;
}

    /*
    min and max Functions
    Simple utility functions to return the smaller or larger of two numbers.
    */
double min(double a, double b){
    if(a<b){return a;}
    return b;
}
double max(double a, double b){
    if(a<b){return b;}
    return a;
}

//Is a 2100Hz Signal present?
void detect2100(vector<complex<double>>* sampleptr){
    vector<complex<double>> samples = *sampleptr;
    double duration=1.0;
    double FRQ=2100/duration;
    int n=samples.size();
    if(n&(n-1)){
        //ignore last values.
        int m=2;
        while (m<n){m=m<<1;}
        m=m>>1;
        for(int i=m;i<n;i++){
            samples.pop_back();
        }
        n=m;
    }
    vector<complex<double>>res=f(samples);
    vector<double> mag(n/2);

    double last=0;
    double llast=0;
    double lllast=0;
    double llllast=0;
    double t=0;
    double floor=0;
    int cx=0;
    for (int i = 1; i < n/2;i++) {
        t=(res[i].real()*res[i].real())+(res[i].imag()*res[i].imag());
        if(i>2){
            mag[i-2]=t+last+llast+lllast+llllast;
        }
        if((i+5>=FRQ && i-5<=FRQ) ){}
        else{
            if(i<1000){
                floor+=t;
                cx++;
            }
        }
        llllast=lllast;
        lllast=llast;
        llast=last;
        last=t;
    }

    double threshold=floor/cx;
    /*
        if res==-32768 -> no signal detected.
        -> ratio ~probability  P(FRQ2)=(res/655.36)%+50%.
        -> if positive: more likely FRQ2, if neg: more likely FRQ1.
    */
    printf("The Noise Floor: %f, for the Frequency 2100 is at %f\n",threshold,mag[FRQ]);
    printf("echoSuppressed before %d-------------------------------------\n", echoSuppressed);
    echoSuppressed=(mag[(int)FRQ]>2*threshold); // cast
    if (echoSuppressed) {
        printf("=======================================CONNECTION ESTABLISHED=======================================");
    }
}

//Look if one of the Frequencies are present; If both are, compute a likelihood; If none are, return -32768.
short detect2FRQ(vector<complex<double>> samples,double threshold,double duration,int FRQ1, int FRQ2){
    FRQ1/=duration;
    FRQ2/=duration;
    int n=samples.size();
    if(n&(n-1)){
        int m=2;
        while (m<n){m=m<<1;}
        m=m>>1;
        for(int i=m;i<n;i++){
            samples.pop_back();
        }
        n=m;
    }
    vector<complex<double>>res=f(samples);
    vector<double> mag(n/2);

    double last=0;
    double llast=0;
    double lllast=0;
    double llllast=0;
    double t=0;
    for (int i = 1; i < n/2;i++) {
        t=(res[i].real()*res[i].real())+(res[i].imag()*res[i].imag());
        if(i>2){
            mag[i-2]=t+last+llast+lllast+llllast;
        }
        llllast=lllast;
        lllast=llast;
        llast=last;
        last=t;
    }

    double floor=0;

    /*
        if res==-32768 -> no signal detected.
        -> ratio ~probability  P(FRQ2)=(res/655.36)%+50%.
        -> if positive: more likely FRQ2, if neg: more likely FRQ1.
    */
    u8 N=(mag[(int)FRQ1]>threshold)+2*(mag[(int)FRQ2]>threshold);
    switch(N){
        case 0:
            return -32768;
        case 1:
            return -32767;
        case 2:
            return 32767;
        case 3:
            return (mag[(int)FRQ2]/(mag[(int)FRQ2]+mag[(int)FRQ1])-0.5)*65536;

    }
    return 0;
}

/*
From portaudio library
*/
static int recordCallback( const void *inputBuffer, void *outputBuffer,
        unsigned long framesPerBuffer,
        const PaStreamCallbackTimeInfo* timeInfo,
        PaStreamCallbackFlags statusFlags,
        void *userData ){

    // data is a paTestData datatype
    paTestData *data = (paTestData*)userData;
    const SAMPLE *rptr = (const SAMPLE*)inputBuffer;
    SAMPLE *wptr = &data->recordedSamples[data->frameIndex * NUM_CHANNELS];
    long framesToCalc;
    long i;
    int finished;
    unsigned long framesLeft = data->maxFrameIndex - data->frameIndex;

    (void) outputBuffer; /* Prevent unused variable warnings. */
    (void) timeInfo;
    (void) statusFlags;
    (void) userData;

    if( framesLeft < framesPerBuffer ){
        framesToCalc = framesLeft;
        finished = paComplete;
    }
    else{
        framesToCalc = framesPerBuffer;
        finished = paContinue;
    }

    if( inputBuffer == NULL ){
        for( i=0; i<framesToCalc; i++ ){
            *wptr++ = SAMPLE_SILENCE;  /* left */
            if( NUM_CHANNELS == 2 ) *wptr++ = SAMPLE_SILENCE;  /* right */
        }
    }
    else{
        for( i=0; i<framesToCalc; i++ ){
            *wptr++ = *rptr++;  /* left */
            if( NUM_CHANNELS == 2 ) *wptr++ = *rptr++;  /* right */
        }
    }
    data->frameIndex += framesToCalc;
    return finished;
}

// From portaudio library
/* This routine will be called by the PortAudio engine when audio is needed.
** It may be called at interrupt level on some machines so don't do anything
** that could mess up the system like calling malloc() or free().
*/
static int playCallback( const void *inputBuffer, void *outputBuffer,
      unsigned long framesPerBuffer,
      const PaStreamCallbackTimeInfo* timeInfo,
      PaStreamCallbackFlags statusFlags,
      void *userData ) {

    paTestData *data = (paTestData*)userData;
    SAMPLE *rptr = &data->recordedSamples[data->frameIndex * NUM_CHANNELS];
    SAMPLE *wptr = (SAMPLE*)outputBuffer;
    unsigned int i;
    int finished;
    unsigned int framesLeft = data->maxFrameIndex - data->frameIndex;

    (void) inputBuffer; /* Prevent unused variable warnings. */
    (void) timeInfo;
    (void) statusFlags;
    (void) userData;

    if( framesLeft < framesPerBuffer ){
        // final buffer...
        for( i=0; i<framesLeft; i++ ){
            *wptr++ = *rptr++;  // left
            if( NUM_CHANNELS == 2 ) *wptr++ = *rptr++;//right
        }
        for( ; i<framesPerBuffer; i++ ){
            *wptr++ = 0;  // left
            if( NUM_CHANNELS == 2 ) *wptr++ = 0;  // right
        }
        data->frameIndex += framesLeft;
        finished = paComplete;
    }
    else{
        for( i=0; i<framesPerBuffer; i++ ){
            *wptr++ = *rptr++;  //left
            if( NUM_CHANNELS == 2 ) *wptr++ = *rptr++;  // right
        }
        data->frameIndex += framesPerBuffer;
        finished = paContinue;
    }
    return finished;
}

/*
This function records audio input using the PortAudio library returns a vector of std::complex<double>

Declarations:
Various variables such as PaStreamParameters, PaStream, and paTestData are declared for handling PortAudio functionalities and storing audio data.

Memory Allocation:
Memory is allocated for recordedSamples to hold the audio samples.
Each sample is initialized to 0.

Recording Audio
PortAudio Initialization:
Pa_Initialize() initializes the PortAudio library.
The default input device is selected using Pa_GetDefaultInputDevice().

Stream Configuration:
Input stream parameters are set, including channel count, sample format, and latency.
The stream is opened using Pa_OpenStream() with the recordCallback function to capture audio data.

Recording Process:
The stream is started with Pa_StartStream(), and the program enters a loop to wait while audio data is recorded.
Progress is logged via printf.

Post-Processing
Amplitude Analysis:
The maximum amplitude and average amplitude of the recorded data are calculated. This provides insights into the loudness of the recorded audio.
Audio data is pushed into inputVector after casting each sample to double.

Playback:
The recorded data is optionally played back for debugging purposes using the default output device.

Cleanup:
PortAudio is terminated with Pa_Terminate().
Allocated memory for recordedSamples is freed to prevent memory leaks.
*/
// Adapted code from portaudio library
vector<complex<double>> listen(double NUM_SECONDS=1.0) {
    vector<complex<double>> inputVector;
    PaStreamParameters  inputParameters,
    outputParameters;
    PaStream*           stream;
    PaError             err = paNoError;
    paTestData          data;
    int                 i;
    int                 totalFrames;
    int                 numSamples;
    int                 numBytes;
    SAMPLE              max, val;
    double              average;

    printf("patest_record.c\n"); fflush(stdout);

    data.maxFrameIndex = totalFrames = NUM_SECONDS * SAMPLE_RATE; /* Record for a few seconds. */
    data.frameIndex = 0;
    numSamples = totalFrames * NUM_CHANNELS;
    numBytes = numSamples * sizeof(SAMPLE);
    // INFO these are our data that are recorded
    data.recordedSamples = (SAMPLE *) malloc( numBytes ); /* From now on, recordedSamples is initialised. */
    if( data.recordedSamples == NULL ){
        printf("Could not allocate record array.\n");
        goto done;
    }
    for( i=0; i<numSamples; i++ ){data.recordedSamples[i] = 0;}

    err = Pa_Initialize();
    if( err != paNoError ){
        goto done;
    }

    inputParameters.device = Pa_GetDefaultInputDevice(); /* default input device */
    if (inputParameters.device == paNoDevice) {
        fprintf(stderr,"Error: No default input device.\n");
        goto done;
    }
    inputParameters.channelCount = NUM_CHANNELS;
    inputParameters.sampleFormat = PA_SAMPLE_TYPE;
    inputParameters.suggestedLatency = Pa_GetDeviceInfo( inputParameters.device )->defaultLowInputLatency;
    inputParameters.hostApiSpecificStreamInfo = NULL;

    // Record some audio. --------------------------------------------
    //NULL: &outputParameters.
    //paClipOff: we won't output out of range samples so don't bother clipping them
    err = Pa_OpenStream(&stream,&inputParameters,NULL,SAMPLE_RATE,FRAMES_PER_BUFFER,paClipOff,recordCallback,&data);
    if( err != paNoError){
        goto done;
    }
    err = Pa_StartStream( stream );
    if( err != paNoError ){
        goto done;
    }
    printf("\n===Now recording!!===\n"); fflush(stdout);
    while( (err=Pa_IsStreamActive(stream) ) == 1 ){
        Pa_Sleep(1000);
        printf("index = %d\n", data.frameIndex ); fflush(stdout);
    }
    if( err < 0 ){
        goto done;
    }

    err = Pa_CloseStream( stream );
    if( err != paNoError ){
        goto done;
    }

    // Measure maximum peak amplitude.
    max = 0;
    average = 0.0;

    for( i=0; i<numSamples; i++ ){
        val = data.recordedSamples[i];
        if( val < 0 ){
            val = -val; //ABS
        }
        if( val > max ){
            max = val;
        }
        average += val;
    }

    average = average / (double)numSamples;

    printf( " sample max amplitude = " PRINTF_S_FORMAT " \n ", max );
    printf( " sample average = %lf\n ", average );

    // Write recorded data to a file.
    // This is also used to work with the data for our program
#if WRITE_TO_FILE
{
    FILE  *fid;
    fid = fopen("recorded.raw", "wb");
    if( fid == NULL ){
        printf("Could not open file.");
    }
    else{
        fwrite( data.recordedSamples, NUM_CHANNELS * sizeof(SAMPLE), totalFrames, fid );
        fclose( fid );
        printf("Wrote data to 'recorded.raw'\n");

        // iterate through the stream and push to vector
        for (int i = 0; i < numSamples; i++) {
            inputVector.push_back(data.recordedSamples[i]);
        }
    }
}
#endif

    /* Playback recorded data.  -------------------------------------------- */
    //INFO: this can be left here for DEBUG purposes to check wether it recorded something or not
    data.frameIndex = 0;

    outputParameters.device = Pa_GetDefaultOutputDevice(); /* default output device */
    if (outputParameters.device == paNoDevice) {
        fprintf(stderr,"Error: No default output device.\n");
        goto done;
    }
    outputParameters.channelCount = NUM_CHANNELS;
    outputParameters.sampleFormat =  PA_SAMPLE_TYPE;
    outputParameters.suggestedLatency = Pa_GetDeviceInfo( outputParameters.device )->defaultLowOutputLatency;
    outputParameters.hostApiSpecificStreamInfo = NULL;
done:
    Pa_Terminate();
    if( data.recordedSamples ){       /* Sure it is NULL or valid. */
        free( data.recordedSamples );
    }
    if( err != paNoError ){
        fprintf( stderr, "An error occurred while using the portaudio stream\n" );
        fprintf( stderr, "Error number: %d\n", err );
        fprintf( stderr, "Error message: %s\n", Pa_GetErrorText( err ) );
        err = 1;          /* Always return 0 or 1, but no other return codes. */
        // this is a problem, but it is a standard code from portAudio to return in case of error
    }
    // returns the vector
    return inputVector;
}



/*
This function detects frequencies present in a signal using the FFT output.
Steps

Size Adjustment:
Ensures the input size is a power of 2 by truncating excess samples.

FFT Computation:
Calls the f function to transform the input to the frequency domain.

Magnitude Calculation:
Computes the magnitude of each frequency component.

Thresholding:
Identifies significant frequencies based on a specified threshold.

Frequency Filtering:
Filters out frequencies of interest (e.g., DTMF tones like 697 Hz, 770 Hz).

Debugging Output:
Prints detected frequencies and their amplitudes for analysis.
*/
//duration in seconds
vector<int> detect(vector<complex<double>> samples,double threshold,double duration){
    int n=samples.size();
    if(n&(n-1)){
        //ignore last values.
        int m=2;
        while (m<n){m=m<<1;}
        m=m>>1;
        for(int i=m;i<n;i++){
            samples.pop_back();
        }
        n=m;
    }
    vector<complex<double>>res=f(samples);
    vector<double> mag(n/2);

    double last=0;
    double llast=0;
    double lllast=0;
    double llllast=0;
    double t=0;
    for (int i = 1; i < n/2;i++) {
        t=(res[i].real()*res[i].real())+(res[i].imag()*res[i].imag());
        if(i>2){
            mag[i-2]=t+last+llast+lllast+llllast;
        }
        llllast=lllast;
        lllast=llast;
        llast=last;
        last=t;
    }

    vector<int>outvals;
    double floor=0;

    for(int i=1;i<n/2;i++){
        if (mag[i]>threshold){outvals.push_back(i);}
        else{floor=max(floor,mag[i]);}
    }
    u64 ecx=0;

    vector<int> actual;
    for(int i=0;i<outvals.size();i++){
        //printf("Found %4d Hz\t",outvals[i]);
        if(outvals[i]==697){ecx++;actual.push_back(outvals[i]);}
        if(outvals[i]==770){ecx++;actual.push_back(outvals[i]);}
        if(outvals[i]==852){ecx++;actual.push_back(outvals[i]);}
        if(outvals[i]==941){ecx++;actual.push_back(outvals[i]);}
        if(outvals[i]==1209){ecx++;actual.push_back(outvals[i]);}
        if(outvals[i]==1336){ecx++;actual.push_back(outvals[i]);}
        if(outvals[i]==1477){ecx++;actual.push_back(outvals[i]);}
        if(outvals[i]==1633){ecx++;actual.push_back(outvals[i]);}
        //if (i%10==9){printf("\n");}
    }
    printf("\nFound:%zu,[%llu] Floor:%f, 697Hz @ %f\n", outvals.size(), ecx, floor, mag[697]);
    printf("\nFound:%zu,[%llu] Floor:%f, 770Hz @ %f\n", outvals.size(), ecx, floor, mag[770]);
    printf("\nFound:%zu,[%llu] Floor:%f, 852Hz @ %f\n", outvals.size(), ecx, floor, mag[852]);
    printf("\nFound:%zu,[%llu] Floor:%f, 941Hz @ %f\n", outvals.size(), ecx, floor, mag[941]);
    printf("\nFound:%zu,[%llu] Floor:%f, 1209Hz @ %f\n", outvals.size(), ecx, floor, mag[1209]);
    printf("\nFound:%zu,[%llu] Floor:%f, 1336Hz @ %f\n", outvals.size(), ecx, floor, mag[1336]);
    printf("\nFound:%zu,[%llu] Floor:%f, 1477Hz @ %f\n", outvals.size(), ecx, floor, mag[1477]);
    printf("\nFound:%zu,[%llu] Floor:%f, 1633Hz @ %f\n", outvals.size(), ecx, floor, mag[1633]);


    printf("=> Actual Frequencies:");
    for(int i=0;i<actual.size();i++){
        printf("%dHz\t",actual[i]);
    }
    printf("\n");
    return outvals;
}




u64 ctr=0;
u64 NFALSE=0;
//Return the found DTMF tone assuming one exists.
char DTMF_detect_frame(vector<complex<double>> samples){
    printf("\n");
    int n=samples.size();
    vector<complex<double>> res=f(samples);
    vector<double> mag(8);
    for(int ii=0;ii<8;ii++){
        int i=FREQS[ii];
        mag[ii]=(res[i].real()*res[i].real())+(res[i].imag()*res[i].imag());
    }
    u8 lower=0;
    double vollower=mag[0];
    for(int i=0;i<4;i++){
        if (mag[i]>vollower){lower=i;vollower=mag[i];}
    }
    u8 upper=0;
    double volupper=mag[4];
    for(int i=4;i<7;i++){//IGNORE LAST UPPER (A-D) FOR NOW.
        if (mag[i]>volupper){upper=i-4;volupper=mag[i];}
    }
    char matrix[]="147*2580369#";
    char correct[]="x ";
    return matrix[(upper*4+lower)];
}


//Detect a sequence of numbers using the function above.
vector<char> DTMF_detect(vector<complex<double>> samples, int window){
    vector<char> result;
    int m=samples.size()/window;
    printf("output:\nx y ->N");
    for(int i=0;i<m;i++){
        vector<complex<double>> frame;
        for(int j=window*i;j<window*i+window;j++){
            frame.push_back(samples[j]);
        }
        char t=DTMF_detect_frame(frame);
        result.push_back(t);
    }
    printf("\n");
    return result;
}

vector<int> detectionOfFrequencies(void){

    vector<complex<double>> dataVector;
    vector<int> outPutVector;
    vector<complex<double>> dataVector2;
    dataVector = listen();

    for (int i = 0; i < 262144; i++) {
        dataVector2.push_back(dataVector[i]);
    }

    int sizeOut = outPutVector.size();
    printf("Size of vector output ");
    cout << sizeOut << endl;

    int sizeIn = dataVector.size();
    printf("Size of vector input ");
    cout << sizeIn << endl;

    cout << "--------------------second detect--------------------" << endl;
    outPutVector = detect(dataVector2, 20.0, 6);

    sizeOut = outPutVector.size();
    printf("Size of vector output ");
    cout << sizeOut << endl;

    sizeIn = dataVector.size();
    printf("Size of vector input ");
    cout << sizeIn << endl;

    return outPutVector;
}

/*-------------------------------END DETCTION SECTION HERE-----------------------------------*/

/*-------------------------------START TRANSMITTER SECTION HERE-----------------------------------*/
// Self tailored Audio callback function for two frequencies from our frequency pool
static int audioCallback2(const void *inputBuffer, void *outputBuffer,
        unsigned long framesPerBuffer,
        const PaStreamCallbackTimeInfo *timeInfo,
        PaStreamCallbackFlags statusFlags,
        void *userData) {

    float *out = (float *)outputBuffer;
    static float phase1 = 0.0f;
    static float phase2 = 0.0f;

    // Retrieve frequencies from userData
    FrequencyPair *frequencies = (FrequencyPair *)userData;

    for (unsigned long i = 0; i < framesPerBuffer; i++) {
        float sine1 = sinf(phase1);
        float sine2 = sinf(phase2);
        float sample = sine1 + sine2;

        *out++ = sample; // Left channel
        *out++ = sample; // Right channel
        phase1 += 2.0f * M_PI * frequencies->frequency1 / SAMPLE_RATE;
        phase2 += 2.0f * M_PI * frequencies->frequency2 / SAMPLE_RATE;

        if (phase1 >= 2.0f * M_PI) phase1 -= 2.0f * M_PI;
        if (phase2 >= 2.0f * M_PI) phase2 -= 2.0f * M_PI;
    }

    return paContinue;
}

/*
Echo suppression output: double vector w/ samples.
(Duration: 1.5s@ 2100Hz, 180deg PhR every 0.45s)
*/
vector<double> echoSuppressionSuppression() {
    vector<double> echosuppressionOut;
    int Ntot = SAMPLERATE * 1.5;  // Total samples for 1.5 seconds
    double dt = 1.0 / SAMPLE_RATE;  // Time step
    double amplitude = 0.5;  // Normalized amplitude to prevent clipping
    for (int i = 0; i < Ntot; i++) {
        double phaseShift = (-1 + 2 * ((i / (int)(0.45 * SAMPLE_RATE)) % 2)); // Alternates every 0.45s
        double sample = amplitude * phaseShift * sin(2.0 * M_PI * 2100 * i * dt);
        echosuppressionOut.push_back(sample);
    }

    printf("Echo Suppression Sound Generated: %lu samples\n", echosuppressionOut.size());

    return echosuppressionOut;
}


/* Tailored audiocallbackFunction for echosuppresion because it was necessary to generate different tones compared to our
frequency matrix we already created for audiocallback2()
*/
static int audioCallbackWithoutFrequencies(const void *inputBuffer, void *outputBuffer,
    unsigned long framesPerBuffer,
    const PaStreamCallbackTimeInfo *timeInfo,
    PaStreamCallbackFlags statusFlags,
    void *userData) {

    float *out = (float *)outputBuffer;
    static unsigned long sampleIndex = 0;

    vector<double> *sampleVectorEcho = (vector<double> *)userData; // Retrieve data

    for (unsigned long i = 0; i < framesPerBuffer; i++) {
        if (sampleIndex >= sampleVectorEcho->size()) {
            sampleIndex = 0; // Restart sound when we reach the end
        }

        float sample = (*sampleVectorEcho)[sampleIndex++];
        *out++ = sample; // Left channel
        *out++ = sample; // Right channel
    }
    return paContinue;
}


/* Function to play a beep for given frequencies for the echosuppression sound.
Separate beeps() function for the setup with the echosuppressionsuppression()
*/
static void beepsSetup() {
    vector<double> bufferedData = echoSuppressionSuppression();

    PaError err = Pa_Initialize();
    if (err != paNoError) {
        cout << "PortAudio error: " << Pa_GetErrorText(err) << endl;
        return;
    }

    PaStream *stream;
    err = Pa_OpenDefaultStream(&stream,0,NUM_CHANNELS,paFloat32,SAMPLE_RATE,FRAMES_PER_BUFFER,audioCallbackWithoutFrequencies,&bufferedData);
    if (err != paNoError) {
        cout << "PortAudio error: " << Pa_GetErrorText(err) << endl;
        return;
    }

    err = Pa_StartStream(stream);
    if (err != paNoError) {
        cout << "PortAudio error: " << Pa_GetErrorText(err) << endl;
        return;
    }

    Pa_Sleep(1500); // Play sound for 1.5 seconds
    // This is always important for durations of sounds to be played

    err = Pa_StopStream(stream);
    if (err != paNoError) {
        cout << "PortAudio error: " << Pa_GetErrorText(err) << endl;
    }

    Pa_CloseStream(stream);
    Pa_Terminate();
}



/* Classic beep for our dial numbers
*/
static void beeps(FrequencyPair frequencies) {
    PaError err = Pa_Initialize();
    if (err != paNoError) {
        cout << "PortAudio error: " << Pa_GetErrorText(err) << endl;
        return;
    }

    PaStream *stream;
    err = Pa_OpenDefaultStream(&stream,0,NUM_CHANNELS,paFloat32,SAMPLE_RATE,FRAMES_PER_BUFFER,audioCallback2,&frequencies);
    if (err != paNoError) {
        cout << "PortAudio error: " << Pa_GetErrorText(err) << endl;
        return;
    }

    err = Pa_StartStream(stream);
    if (err != paNoError) {
        cout << "PortAudio error: " << Pa_GetErrorText(err) << endl;
        return;
    }

    Pa_Sleep(DURATION * 1000); // Play sound for the duration
    // This is always important for durations of sounds to be played

    err = Pa_StopStream(stream);
    if (err != paNoError) {
        cout << "PortAudio error: " << Pa_GetErrorText(err) << endl;
    }

    Pa_CloseStream(stream);
    Pa_Terminate();
}

/*Function to play a sequence of tones based on input numbers from the frequency matrix
*/
static void beepSequenceOfNumbers(vector<u8> inputNumberSequence, FrequencyPair frequencies) {
    for (u8 number : inputNumberSequence) {
        if (number >= '0' && number <= '9') {
            int digit = number - '0'; // Convert char to digit
            frequencies.frequency1 = frqU[digit];
            frequencies.frequency2 = frqD[digit];
            beeps(frequencies);
        } else {
            cout << "Invalid input detected: " << number << endl;
        }
    }
}

/*-------------------------------START V.8. BIS SECTION HERE-----------------------------------*/

// Function to generate a preamble for the communication sequence
vector<double> preamble(){
    float dt = 1.0f / SAMPLERATE;  // Time step based on sample rate
    float dsSamples = SAMPLERATE / 10;  // Duration of preamble in samples
    vector<double> out;

    // Generate sine wave for the preamble signal at 980Hz
    for(int i = 0; i < dsSamples; i++){
        out.push_back(sin(980 * dt * i * THETA));  // Sine wave generation
    }

    return out;
}


//encode data using the V.21 standard.
//isremote = 0/1, 0: the calling one, 1: the receiving one
//number: input to send
//data: samples which get sent, aka output

void encodeData(vector<double>* data, u8 number, u8 isremote){
    int Nsamples = SAMPLERATE / 30;  // Duration for each bit at 30bps // V.21 standard is 300bps
    double dt = 1.0f / SAMPLERATE;  // Time step, 1 sample has this much time in it

    // Process each bit (8 bits in total)
    for(int ii = 7; ii >= 0; ii--){
        u8 i = (number >> ii) & 1;  // Get each bit
        float FRQ = 980 + 200 * i + 670 * isremote;  // Choose frequency based on bit value and remote/local flag. Defines who's sending
        for(int j = 0; j < Nsamples; j++){
            data->push_back(sin(FRQ * j * dt * THETA));  // Encode the bit as a sine wave at the chosen frequency
        }
    }
}

// Function to initiate a call (sent from the remote side)
vector<double> initiateCallA(){
    float dt = 1.0f / SAMPLERATE;
    float dsSamples = SAMPLERATE / 10;
    vector<double> out = preamble();  // Start with the preamble

    // Generate MSe signal (Mixed signal at 1375Hz and 2002Hz)
    // The tones come one after the other, to announce the willing of sending MSe signal from the v8bis standard:
    // MSe: Mode select Send, MSe from local to receiving end. Call establishement
    // one tone is specific
    for(int i = 0; i < dsSamples * 4; i++){
        out.push_back(sin(1375 * dt * i * THETA) * 0.5 + sin(2002 * dt * i * THETA) * 0.5);
    }

    // Generate 650Hz signal
    for(int i = 0; i < dsSamples * 1; i++){
        out.push_back(sin(650 * dt * i * THETA));
    }

    // Send ModeRequest (5 flags: 01111110)
    printf("===> ENCODEDATA: LENGTH=%d",out.size());
    for(int i = 0; i < 5; i++) {
        // Flags to synchronise
        encodeData(&out, 0b01111110, 1);  // Mode Request Flag
        //changed to dereferenced because it throws an error in the intellisense analysis of vscode
    }
    printf("ENCODEDATA =>AFTER FLAGS=%d",out.size());

    // Npar = 0, Spar1 = 1.0, Npar2 = 0(1), Spar2 = 0.8.0

    // Send specific data (MS, parameters, checksum)
    // changed to dereferenced because it throws an error in the intellisense analysis of vscode
    encodeData(&out, 0b00100001, 1);  // Mode Set (MS): Mode select

    encodeData(&out, 0b00000000, 1);  // Npar1: Parameter 1 with single value. Here No parameters = 0

    encodeData(&out, 0b00000001, 1);  // Spar1a: Sub-Parameter 1 with mutliple parts to it. Here Network type = 1, meaning Data only
    encodeData(&out, 0b00000000, 1);  // Spar1b: Adjacent to a.)

    encodeData(&out, 0b00000000, 1);  // Npar2a: Analogous to Spar. Here no parameter = 0.
    encodeData(&out, 0b00000000, 1);  // Npar2b: Here meaning file transfer no parameter = 0
    encodeData(&out, 0b00001000, 1);  // Npar2c: Speed mode V.21

    // Checksum (0x7ffb): polynomial: (x^15 - x^14 - ... x^0)
    // x = Npar1 Spar1a Spar1b Npar2a Npar2b Npar2C: stuck together as a binary
    
    encodeData(&out, 0x7f, 1);
    encodeData(&out, 0xfb, 1);


    // Send final flags (end of sequence)
    // Send flags as above
    for(int i = 0; i < 3; i++) {
        encodeData(&out, 0b01111110, 1);  // Mode Request Flag
    }

    printf("Initiating V.21 Communication via V.8bis\n");
    return out;
}

// Function to receive a call (sent from the local side)
vector<double> receiveCallA(){
    vector<double> out = preamble();  // Start with the preamble
    float dt = 1.0f / SAMPLERATE;
    float dsSamples = SAMPLERATE / 10;

    // Generate ESCAPE signal (1529Hz + 2225Hz): Escape sequence is for the defining the behaviour to sent messages
    for(int i = 0; i < dsSamples * 4; i++){
        out.push_back(sin(1529 * dt * i * THETA) * 0.5 + sin(2225 * dt * i * THETA) * 0.5);
    }

    // Generate 1650Hz signal
    for(int i = 0; i < dsSamples * 1; i++){
        out.push_back(sin(1650 * dt * i * THETA));
    }

    // Send ModeRequest (5 flags: 01111110)
    for(int i = 0; i < 5; i++) {
        encodeData(&out, 0b01111110, 0);  // Mode Request Flag
    }

    // Send specific data (MS, parameters, checksum)
    // Analogous to above
    encodeData(&out, 0b00100001, 0);  // Mode Set (MS)

    encodeData(&out, 0b00000000, 0);  // Npar1

    encodeData(&out, 0b00000001, 0);  // Spar1a
    encodeData(&out, 0b00000000, 0);  // Spar1b

    encodeData(&out, 0b00000000, 0);  // Npar2a
    encodeData(&out, 0b00000000, 0);  // Npar2b
    encodeData(&out, 0b00001000, 0);  // Npar2c


    // Checksum (0x0014)
    encodeData(&out, 0x00, 0);
    encodeData(&out, 0x14, 0);


    // Send final flags (end of sequence)
    for(int i = 0; i < 3; i++) {
        encodeData(&out, 0b01111110, 0);  // Mode Request Flag
    }

    printf("Awaiting Acknowledgement for V.21 communication\n");
    return out;
}

// Function to initiate a call (sent from the local side, ACK)
// Answer to receiveA()
vector<double> initiateCallB(){
    vector<double> out = preamble();  // Start with the preamble

    // Send 5 flags to start the sequence
    for(int i = 0; i < 5; i++) {
        encodeData(&out, 0b01111110, 1);  // Mode Request Flag
    }

    // Send ACK
    encodeData(&out, 0b00100100, 1);  // ACK(1)

    encodeData(&out, 0b00000000, 1);  // Parameters = 0
    encodeData(&out, 0b00000000, 1);
    encodeData(&out, 0b00000000, 1);

    // Checksum (0xf3af)
    encodeData(&out, 0xf3, 1);
    encodeData(&out, 0xaf, 1);

    // Send final flags (end of sequence)
    for(int i = 0; i < 3; i++) {
        encodeData(&out, 0b01111110, 1);  // Mode Request Flag
    }

    printf("Ready for V.21 Communication\n");
    return out;
}

u8 detectPreamble(vector<complex<double>>* sampleptr){
    vector<complex<double>> samples = *sampleptr;
    double duration=0.1;
    double FRQ=980/duration;
    int n=samples.size();
    if(n&(n-1)){
        //ignore last values.
        int m=2;
        while (m<n){m=m<<1;}
        m=m>>1;
        for(int i=m;i<n;i++){
            samples.pop_back();
        }
        n=m;
    }
    vector<complex<double>>res=f(samples);
    vector<double> mag(n/2);

    double last=0;
    double llast=0;
    double lllast=0;
    double llllast=0;
    double t=0;
    double floor=0;
    int cx=0;
    for (int i = 1; i < n/2;i++) {
        t=(res[i].real()*res[i].real())+(res[i].imag()*res[i].imag());
        if(i>2){
            mag[i-2]=t+last+llast+lllast+llllast;
        }
        if((i+5>=FRQ && i-5<=FRQ)){}
        else{
            if(i<1000){
                floor+=t;
                cx++;
            }
        }
        llllast=lllast;
        lllast=llast;
        llast=last;
        last=t;
    }

    double threshold=floor/cx;
    return ((mag[(int)FRQ]>threshold));

}

u8 detectMSe1(vector<complex<double>>* sampleptr){
    vector<complex<double>> samples = *sampleptr;
    double duration=0.4;
    double FRQ1=1375/duration;
    double FRQ2=2002/duration;
    int n=samples.size();
    if(n&(n-1)){
        //ignore last values.
        int m=2;
        while (m<n){m=m<<1;}
        m=m>>1;
        for(int i=m;i<n;i++){
            samples.pop_back();
        }
        n=m;
    }
    vector<complex<double>>res=f(samples);
    vector<double> mag(n/2);

    double last=0;
    double llast=0;
    double lllast=0;
    double llllast=0;
    double t=0;
    double floor=0;
    int cx=0;
    for (int i = 1; i < n/2;i++) {
        t=(res[i].real()*res[i].real())+(res[i].imag()*res[i].imag());
        if(i>2){
            mag[i-2]=t+last+llast+lllast+llllast;
        }
        if((i+5>=FRQ1 && i-5<=FRQ1) ||(i+5>=FRQ2 && i-5<=FRQ2) ){}
        else{
            if(i<1000){
                floor+=t;
                cx++;
            }
        }
        llllast=lllast;
        lllast=llast;
        llast=last;
        last=t;
    }

    double threshold=floor/cx;
    return ((mag[(int)FRQ1] + mag[(int)FRQ2])>threshold);

}

u8 detectMSe2(vector<complex<double>>* sampleptr){
    vector<complex<double>> samples = *sampleptr;
    double duration=0.1;
    double FRQ=650/duration;
    int n=samples.size();
    if(n&(n-1)){
        //ignore last values.
        int m=2;
        while (m<n){m=m<<1;}
        m=m>>1;
        for(int i=m;i<n;i++){
            samples.pop_back();
        }
        n=m;
    }
    vector<complex<double>>res=f(samples);
    vector<double> mag(n/2);

    double last=0;
    double llast=0;
    double lllast=0;
    double llllast=0;
    double t=0;
    double floor=0;
    int cx=0;
    for (int i = 1; i < n/2;i++) {
        t=(res[i].real()*res[i].real())+(res[i].imag()*res[i].imag());
        if(i>2){
            mag[i-2]=t+last+llast+lllast+llllast;
        }
        if((i+5>=FRQ && i-5<=FRQ)){}
        else{
            if(i<1000){
                floor+=t;
                cx++;
            }
        }
        llllast=lllast;
        lllast=llast;
        llast=last;
        last=t;
    }

    double threshold=floor/cx;
    return ((mag[(int)FRQ]>threshold));

}


u8 detectESr1(vector<complex<double>>* sampleptr){
    vector<complex<double>> samples = *sampleptr;
    double duration=0.4;
    double FRQ1=1529/duration;
    double FRQ2=2225/duration;
    int n=samples.size();
    if(n&(n-1)){
        //ignore last values.
        int m=2;
        while (m<n){m=m<<1;}
        m=m>>1;
        for(int i=m;i<n;i++){
            samples.pop_back();
        }
        n=m;
    }
    vector<complex<double>>res=f(samples);
    vector<double> mag(n/2);

    double last=0;
    double llast=0;
    double lllast=0;
    double llllast=0;
    double t=0;
    double floor=0;
    int cx=0;
    for (int i = 1; i < n/2;i++) {
        t=(res[i].real()*res[i].real())+(res[i].imag()*res[i].imag());
        if(i>2){
            mag[i-2]=t+last+llast+lllast+llllast;
        }
        if((i+5>=FRQ1 && i-5<=FRQ1) ||(i+5>=FRQ2 && i-5<=FRQ2) ){}
        else{
            if(i<1000){
                floor+=t;
                cx++;
            }
        }
        llllast=lllast;
        lllast=llast;
        llast=last;
        last=t;
    }

    double threshold=floor/cx;
    return ((mag[(int)FRQ1] + mag[(int)FRQ2])>threshold);

}


u8 detectESr2(vector<complex<double>>* sampleptr){
    vector<complex<double>> samples = *sampleptr;
    double duration=0.1;
    double FRQ=1650/duration;
    int n=samples.size();
    if(n&(n-1)){
        //ignore last values.
        int m=2;
        while (m<n){m=m<<1;}
        m=m>>1;
        for(int i=m;i<n;i++){
            samples.pop_back();
        }
        n=m;
    }
    vector<complex<double>>res=f(samples);
    vector<double> mag(n/2);

    double last=0;
    double llast=0;
    double lllast=0;
    double llllast=0;
    double t=0;
    double floor=0;
    int cx=0;
    for (int i = 1; i < n/2;i++) {
        t=(res[i].real()*res[i].real())+(res[i].imag()*res[i].imag());
        if(i>2){
            mag[i-2]=t+last+llast+lllast+llllast;
        }
        if((i+5>=FRQ && i-5<=FRQ)){}
        else{
            if(i<1000){
                floor+=t;
                cx++;
            }
        }
        llllast=lllast;
        lllast=llast;
        llast=last;
        last=t;
    }

    double threshold=floor/cx;
    return ((mag[(int)FRQ]>threshold));

}

u8 decodeA(vector<complex<double>>* buf) {
    vector<complex<double>> fragment;
    int fraglength = SAMPLERATE/30;
    vector<complex<double>> buffer = *buf;
    
    double silenceThresh=20;

    // Find the first significant signal by the manually setting the threshold
    int startIdx;
    for (startIdx = 0; startIdx < buffer.size(); startIdx++) {
        if (abs(buffer[startIdx]) > silenceThresh) {
            break;
        }
    }


    // Remove initial silent part
    buffer = vector<complex<double>>(buffer.begin() + startIdx, buffer.end());

    vector<short> resultLikelihoods;
        double acc=0;
        for(int j=0;j<buffer.size();j++){
                acc+=buffer[j].real(); // make it a real number
        }
        double threshold=acc/buffer.size();
        short tmp;
        int cx=0;

        fragment.clear();
        if(cx+226380>fragment.size()){return 0;}
        for(int i=0;i<SAMPLE_RATE/10;i++){
            fragment.push_back(buffer[cx++]);
        }
        if(detectPreamble(&fragment)!=1){return 0;}
        fragment.clear();

        for(int i=0;i<(4*SAMPLE_RATE)/10;i++){
            fragment.push_back(buffer[cx++]);
            
        }
        
        if(detectMSe1(&fragment)!=1){return 0;}
        fragment.clear();
        for(int i=0;i<SAMPLE_RATE/10;i++){
            fragment.push_back(buffer[cx++]);
        }
        if(detectMSe2(&fragment)!=1){return 0;}
        fragment.clear();

        for(int i=0;i<17*8;i++){
            fragment.clear();
            for(int j=0;j<fraglength;j++){
                fragment.push_back(buffer[cx++]);
            }
            tmp=detect2FRQ(fragment,threshold,2,V21FRQ[issender*2+0],V21FRQ[issender*2+1]);
            resultLikelihoods.push_back(tmp);
            printf("%d ",tmp);
        }
    
    vector<char> bts;
    for(int i=0;i<17*8;i++){
        bts.push_back(resultLikelihoods[i]>0);
    }
    unsigned char wanted[]={0x7e,0x7e,0x7e,0x7e,0x7e,
        0x32,0x00,2,0,0,
        0,8,127,0xfb,
        0x7e,0x7e,0x7e};
    int dx=0;
    u8 tmp2=0;
    for(int i=0;i<17*8;i++){
        if(i%8==0){
            tmp2=0;
            if(i){
                if(tmp2!=wanted[dx]){return 0;}
            }
        }
        tmp2*=2;
        tmp2+=bts[i];
    }
    return 1;

}

u8 decodeB(vector<complex<double>>* buf) {
    vector<complex<double>> fragment;
    int fraglength = SAMPLERATE/30;
    vector<complex<double>> buffer = *buf;
    double silenceThresh=20;

    // Find the first significant signal by the manually set threshold
    int startIdx;
    for (startIdx = 0; startIdx < buffer.size(); startIdx++) {
        if (abs(buffer[startIdx]) > silenceThresh) {
            break;
        }
    }


    // Remove initial silent part
    buffer = vector<complex<double>>(buffer.begin() + startIdx, buffer.end());

    vector<short> resultLikelihoods;
        double acc=0;
        for(int j=0;j<buffer.size();j++){
                acc+=buffer[j].real(); // make it a real number
        }
        double threshold=acc/buffer.size();
        short tmp;
        int cx=0;
        if(cx+226380>fragment.size()){return 0;}
        fragment.clear();
        for(int i=0;i<SAMPLE_RATE/10;i++){
            fragment.push_back(buffer[cx++]);
        }
        if(detectPreamble(&fragment)!=1){return 0;}

        if(cx+169050>fragment.size()){return 0;}
        fragment.clear();

        for(int i=0;i<14*8;i++){
            fragment.clear();
            for(int j=0;j<fraglength;j++){
                fragment.push_back(buffer[cx++]);
            }
            tmp=detect2FRQ(fragment,threshold,2,V21FRQ[issender*2+0],V21FRQ[issender*2+1]);
            resultLikelihoods.push_back(tmp);
            printf("%d ",tmp);
        }
    
    vector<char> bts;
    for(int i=0;i<14*8;i++){
        bts.push_back(resultLikelihoods[i]>0);
    }
    unsigned char wanted[]={0x7e,0x7e,0x7e,0x7e,0x7e,
        0x24,0,0,0,
        0xf3,0xaf,
        0x7e,0x7e,0x7e
    };
    int dx=0;
    u8 tmp2=0;
    for(int i=0;i<14*8;i++){
        if(i%8==0){
            tmp2=0;
            if(i){
                if(tmp2!=wanted[dx]){return 0;}
            }
        }
        tmp2*=2;
        tmp2+=bts[i];
    }
    return 1;

}
u8 decodeR(vector<complex<double>>* buf) {
    vector<complex<double>> fragment;
    int fraglength = SAMPLERATE/30;
    vector<complex<double>> buffer = *buf;


    double silenceThresh=20;

    // Find the first significant signal by the manually set threshold
    int startIdx;
    for (startIdx = 0; startIdx < buffer.size(); startIdx++) {
        if (abs(buffer[startIdx]) > silenceThresh) {
            break;
        }
    }

    // Remove initial silent part
    buffer = vector<complex<double>>(buffer.begin() + startIdx, buffer.end());

    vector<short> resultLikelihoods;
        double acc=0;
        for(int j=0;j<buffer.size();j++){
                acc+=buffer[j].real(); // make it a real number
        }
        double threshold=acc/buffer.size();
        short tmp;
        int cx=0;
        if(cx+226380>fragment.size()){return 0;}
        fragment.clear();
        for(int i=0;i<SAMPLE_RATE/10;i++){
            fragment.push_back(buffer[cx++]);
        }
        if(detectPreamble(&fragment)!=1){return 0;}


        for(int i=0;i<(4*SAMPLE_RATE)/10;i++){
            fragment.push_back(buffer[cx++]);
            
        }
        fragment.clear();
        if(detectESr1(&fragment)!=1){return 0;}
        fragment.clear();
        for(int i=0;i<SAMPLE_RATE/10;i++){
            fragment.push_back(buffer[cx++]);
        }
        if(detectESr2(&fragment)!=1){return 0;}
        fragment.clear();

        for(int i=0;i<17*8;i++){
            fragment.clear();
            for(int j=0;j<fraglength;j++){
                fragment.push_back(buffer[cx++]);
            }
            tmp=detect2FRQ(fragment,threshold,2,V21FRQ[issender*2+0],V21FRQ[issender*2+1]);
            resultLikelihoods.push_back(tmp);
            printf("%d ",tmp);
        }
    
    vector<char> bts;
    for(int i=0;i<17*8;i++){
        bts.push_back(resultLikelihoods[i]>0);
    }
    unsigned char wanted[]={0x7e,0x7e,0x7e,0x7e,0x7e,
        0x21,0,1,0,0,0,1,0,0x14,
        0x7e,0x7e,0x7e
    };
    int dx=0;
    u8 tmp2=0;
    for(int i=0;i<17*8;i++){
        if(i%8==0){
            tmp2=0;
            if(i){
                if(tmp2!=wanted[dx]){return 0;}
            }
        }
        tmp2*=2;
        tmp2+=bts[i];
    }
    return 1;

}
void send(vector<double>*);
//communication error: terminate communication.
void crash() {
    vector<double> out=preamble();
    encodeData(&out,0xb1,issender);
    encodeData(&out,0x7f,issender);
    encodeData(&out,0xa7,issender);
#ifdef IGNORE_ERRORS
    send(&out);
    exit(1);
#endif

}
void RX(int);
void setup(int sender){
    cout << "Sender: ";
    cout << sender;
    if(sender){
        // Dial up the number
        // 1. step
        // Hardcoded "phone number" to dial to the receiver part
        cout << "Call telephone number for the receiver: " << endl;
        vector<u8> inputNumberSequence = {'0', '6', '1', '0', '7', '9'};
        cout << "Hardocded number is: 061079" << endl;
        beepSequenceOfNumbers(inputNumberSequence, frequencies);
        // Preamble to give the signal to suppress echo when calling up a conventional modem
        // echosuppression sound dial up generated by the code
        cout << "Echosupression Call: " << endl;
        echosuppressionOut = echoSuppressionSuppression();
        // Separate beeps() function for the setup with the echosuppressionsuppression()
        beepsSetup();

        // Receive call
        vector<complex<double>> buf = listen(9);
        
        if (decodeA(&buf) == 0) {crash();} // Calls a message that it received an error and then exits
        // 0x7e7e7e7e7e 0x3200 0200 0000 087f fb 0x7e7e7e
        
        // after initiatecallA next step here
        // 3. step
        vector<double> bufRecA = receiveCallA();
        send(&bufRecA);
        buf.clear();
        buf = listen(7);
        
        if (decodeB(&buf) == 0) {crash();}; // Calls a message that it received an error and then exits
        // 0x7e7e7e7e7e 0x2400 0000 f3af 0x7e7e7e
        cout << "Setup step for sender done enter sending mode..." << endl;

    }else{
        //detect echosuppression suppression
        RX(1);

        // step 2.
        vector<double> bufInitA = initiateCallA();
        send(&bufInitA);
        vector<complex<double>> buf = listen(7); // These seconds are predefined
        if (decodeR(&buf) == 0) {crash();}; // Calls a message that it received an error and then exits
        // 0x7e7e7e7e7e 0x210001000000080014 0x7e7e7e

        // 4. step
        vector<double> bufB  = initiateCallB();
        send(&bufB);
        cout << "Setup step for sender done enter receiving mode..." << endl;
    }
}


std::atomic<bool> rxThreadBusy(false);
vector<char> finalMessageBits;
void RXdecode(vector<complex<double>>* buf) {
    char key[4]="-01";
    vector<complex<double>> buffer = *buf;
        int fraglength=SAMPLERATE/30;// Length of one symbol (10% of V.21)
        vector<complex<double>> fragment;
        vector<short> resultLikelihoods;
        double acc=0;
        
        for(int j=0;j<buffer.size();j++){
                acc+=buffer[j].real(); // make it a real number
        }
        double threshold=acc/buffer.size();
        short tmp;
        for(int i=0;i<30;i++){
            fragment.clear();
            for(int j=0;j<fraglength;j++){
                fragment.push_back(buffer[i*fraglength+j]);
            }
            tmp=detect2FRQ(fragment,threshold,2,V21FRQ[issender*2+0],V21FRQ[issender*2+1]);
            resultLikelihoods.push_back(tmp);
            printf("%d ",tmp);
        }
        vector<char> txt;
        char t=0;
        for(int i=0;i<resultLikelihoods.size();i++){
            finalMessageBits.push_back(key [(resultLikelihoods[i]!=-32768)*(1+(resultLikelihoods[i]>0))  ]);
            t=(t*2+(resultLikelihoods[i]>0));// ML estimator, best one if Data has high Entropy.
            if(i%8==7){
                txt.push_back(t);
                t=0;    
            }
        }
        printf("The fragment is:");
        for(int i=0;i<txt.size();i++){
            printf("%c",txt[i]);
        }
        printf("\n");
        

        printf("===>Message so far: ");
        for(int i=0;i<finalMessageBits.size();i++){
            printf("%c",finalMessageBits[i]);
            if(i%8==7){printf(" ");}
        }
        printf("\n");
        return;
        }

// Here CHATGPT was used:
// For the thread intern workings I(G.U.) asked chatgpt to help me because I did not work with threads a lot in C++.
// Use of CHATGPT to generate pseudocode based on the concept provided by R.C.
void RX(int isstart) {
    cout << "RX started.\n";
    vector<vector<complex<double>>> RXbuf(2);
    int j = 0;
    future<void> rxFuture;
    while (true) {
        if (echoSuppressed){echoSuppressed=0;return;}
        RXbuf[j].clear();
        cout << "while(1) in RX loop...\n";
        RXbuf[j] = listen(2);  // Capture audio for 2 seconds and put it into the rx buffer we listen for 2sec
        if (echoSuppressed){echoSuppressed=0;return;}
        // Check if the previous async task is still running
        if (rxFuture.valid() && rxFuture.wait_for(chrono::seconds(0)) != future_status::ready) {
            cout << "Previous processing still running, skipping new task.\n";
            continue;  // Skip starting a new thread if the previous one is still running
        }

        // Start a new async task for processing
        if (echoSuppressed){echoSuppressed=0;return;}
        cout << "Starting new processing task...\n";
        if (isstart) {
            if (echoSuppressed){echoSuppressed=0;return;}
            cout << "Calling detect2100()...\n";
            rxFuture = async(launch::async, detect2100, &RXbuf[j]);
            if (echoSuppressed){echoSuppressed=0;return;}
        } else {
            cout << "Calling RXdecode()...\n";
            rxFuture = async(launch::async, RXdecode, &RXbuf[j]);
        }

        j ^= 1;  // Swap buffers
    }
}

std::atomic<bool> txThreadBusy(false);

/* Send is a generalised function to send double vector via loudspeaker, analogous to the beepsSequence and beepsSetup
*/
void send(vector<double>* buf) {
    cout << "-------------------------------------------- START send method -------------------------------------------- " << endl;
    vector<double> bufferedData = *buf;
    cout << "Size of bufferedData vector: ";
    cout << bufferedData.size() << endl;

    double mx=0;
    for(int i=0;i<bufferedData.size();i++){
        mx=max(bufferedData[i]*bufferedData[i],mx);
    }
    mx=sqrt(mx);
    for(int i=0;i<bufferedData.size();i++){
        bufferedData[i]/=mx;
    }

    int durationOfSend = (bufferedData.size()*1000) / SAMPLE_RATE; // 1000 because of milliseconds

    PaError err = Pa_Initialize();
    if (err != paNoError) {
        cout << "PortAudio error: " << Pa_GetErrorText(err) << endl;
        return;
    }

    PaStream *stream;
    err = Pa_OpenDefaultStream(&stream,0,NUM_CHANNELS,paFloat32,SAMPLE_RATE,FRAMES_PER_BUFFER,audioCallbackWithoutFrequencies,&bufferedData);
    if (err != paNoError) {
        cout << "PortAudio error: " << Pa_GetErrorText(err) << endl;
        return;
    }

    err = Pa_StartStream(stream);
    if (err != paNoError) {
        cout << "PortAudio error: " << Pa_GetErrorText(err) << endl;
        return;
    }

    Pa_Sleep(durationOfSend); // Play sound for calculated seconds

    err = Pa_StopStream(stream);
    if (err != paNoError) {
        cout << "PortAudio error: " << Pa_GetErrorText(err) << endl;
    }

    Pa_CloseStream(stream);
    Pa_Terminate();
    cout << "-------------------------------------------- END send method -------------------------------------------- " << endl;
}

void TX(){
    const char txt[]="Hello, I am the system administrator, my voice is my passport, verify me.";
    vector<double> TXbuf;
    for(int i=0;i<74;i++){
        encodeData(&TXbuf,txt[i],0);
    }
    send(&TXbuf);
    printf("%s",txt);
    printf("\n");
}

// Here CHATGPT was used:
// For the thread intern workings I(G.U.) asked chatgpt to help me because I did not work with threads a lot in C++.
// Use of CHATGPT to generate pseudocode based on the concept provided by R.C.
void TX(char txt[], int length) {
    vector<double> TXbuf[2] = {vector<double>(), vector<double>()};
    int Nblocks = length / 40;
    int Nsamples = 1 + (SAMPLE_RATE * 320) / 300; // Adjusted sample calculation
    std::thread txThread;
    int senderFlag = issender;
    int j = 0;
    for (int i = 0; i < Nsamples; i++) {
        TXbuf[0].push_back(0);
        TXbuf[1].push_back(0);
    }
    for (; j < Nblocks - 1; j++) {
        TXbuf[j % 2].clear();
        for (int i = 0; i < 40; i++) {
            encodeData(&TXbuf[j % 2], txt[i + 40 * j], senderFlag);
        }
        if (j && txThread.joinable()) {
            txThread.join(); // Ensure last transmission is done
        }
        // Start new thread for sending
        txThread = std::thread(send, &TXbuf[j % 2]);
    }
    // Process last block
    for (int i = 0; (i + j * 40) < length; i++) {
        encodeData(&TXbuf[j % 2], txt[i + 40 * j], senderFlag);
    }
    if (j && txThread.joinable()) {
        txThread.join();
    }
    txThread = std::thread(send, &TXbuf[j % 2]);
    if (txThread.joinable()) {
        txThread.join(); // Wait for the last transmission to finish
    }
}


/*Tests to test the single functionalities
*/
void testEchoSuppressionOutPutSound() {
    cout << "Test testEchoSuppressionOutPutSound" << endl;
    setup(1);
    cout << "End of Test" << endl;
    return;

}

void testSequenceOfNumbers() {
    const char* inputString = "234";
    size_t length = strlen(inputString);
    vector<u8> inputs;
    for (size_t i = 0; i < length; i++) {
        if (isdigit(inputString[i])) {
            inputs.push_back(inputString[i]);
        } else {
            cout << "Invalid character found: " << inputString[i] << endl;
        }
    }

    FrequencyPair frequencies = frequencies;
    cout << "BEEP SEQUENCE OF NUMBERS" << endl;
    beepSequenceOfNumbers(inputs, frequencies);

    cout << "Done" << endl;

}

void testSendInitA() {
    vector<double> bufA  = initiateCallA();
    send(&bufA);
}

void testSendInitB() {
    vector<double> bufB  = initiateCallB();
    send(&bufB);
}

void testRecCallA() {
    vector<double> bufRecA = receiveCallA();
    send(&bufRecA);
}

// g++ -o main transceiver.cpp -lportaudio -lm
int main() {
    // The test are kind of a status report where we are right now
    /*
    cout << "Test Sequence of numbers in main: " <<endl;
    testSequenceOfNumbers();

    cout << "Test Echossuppression in main: " <<endl;
    testEchoSuppressionOutPutSound();

    cout << "Test initcallA in main: " <<endl;
    testSendInitA();

    cout << "Test initcallB in main: " <<endl;
    testSendInitB();

    cout << "Test receiveCallB in main: " <<endl;
    testRecCallA();
    */

    // main
    printf("Enter 1 to send or 0 to receive: ");
    scanf("%d", &issender);
    issender &= 1;
    cout << "Setup for the transmitter" << endl;
    setup(issender);
    
    // Wait for a few seconds before communicating
    this_thread::sleep_for(chrono::seconds(4));
    if(issender){TX();} else {
        RX(0);}
    return 1;
}
