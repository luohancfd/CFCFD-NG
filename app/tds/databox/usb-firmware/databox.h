/* databox.h */

#ifndef __DATABOX_H
#define __DATABOX_H

// A/D Card Control
void writeCardControlRegister(void);
void resetCardPointer(void);
void select_channel(void);
void select_card(void);
void send_encoded_word(unsigned int);
unsigned int rotate_right(unsigned int c);
void getData(void);
void octalOutput(void);
void voltageOutput(void);
void fullScaleRange(void);
void dataCoupling(void);

// A/D Card Status
unsigned char readCardStatusRegisterHigh(unsigned char n);
unsigned char readCardStatusRegisterLow(unsigned char n);
unsigned char card_is_present(unsigned char n);
unsigned char card_is_sampling(unsigned char n);
unsigned char any_cards_sampling(void);
void report_sampling_status(unsigned char n);
unsigned char whichTimeBase(void);
int getBufferPointer(void);
void printBufferPointer(void);

// Trigger Control
void writeTriggerUnitControlRegister(unsigned char ISAData);
void arm_box(void);
void trigger_box(void);

// Time Base Status
unsigned char readTimeStatusRegisterHigh(unsigned char n);
unsigned char readTimeStatusRegisterLow(unsigned char n);
unsigned char box_is_present(void);
void trigger_number(unsigned char n);
void preTriggerDelay(unsigned char n);
void bufferSizeSelect(unsigned char n);
unsigned char whichTriggerUnit(unsigned char n);
void samplePeriod(unsigned char n);
void timeBaseMultiplier(unsigned char n);

// Trigger Unit Status 
unsigned char readTriggerStatusRegisterHigh(void);
unsigned char readTriggerStatusRegisterLow(void);
void report_coupling(unsigned char n);
void report_slope(unsigned char n);
void reportSlopeAndCouple(unsigned char n);
void enable_threshold_report(unsigned char n);
void report_threshold(unsigned char n);
void report_threshold_in_BCD(unsigned char n);
unsigned char thresholdStatus(void);

// System Functions
void reportIOBase(void);
void setIOBase(void);
void setVerbose(void);
void menu(void);

#endif
