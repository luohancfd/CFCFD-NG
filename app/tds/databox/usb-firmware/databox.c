/* databox.c
 *
 * An set of instructions used specifically for the Databox Interface Card.
 * This file implements the high level user functionality required by the
 * end user.
 * 
 * Luke Hillyard, 15/10/2006
 * Peter Jacobs,  12-Jan-2007
 */

#include <avr/io.h>
#include <avr/eeprom.h>
#include <avr/pgmspace.h>
#include "bits.h"
#include "isa-bus.h"
#include "databox.h"
#include "uart-interrupt.h"

// Addressing Databox Registers
int EEMEM IOBaseAddress   = 0x320;    // 110
int EEMEM TimeBaseAddress = 0x310;

// Most of the operations on the A/D cards work on the 
// currently selected card.
unsigned char CardNum = 1;    // 1 to 15 inclusive
unsigned char ChannelNum = 1; // 1 to 3  inclusive
unsigned char CardReset = 0;
extern unsigned char verbose;

unsigned char suppress_k_report = 1;

// ------------------ A/D Card Control ------------------------

void writeCardControlRegister(void) {
  // Writes byte to the card control register and then sends
  // corresponding text to the serial port.
  unsigned char ISAData = (CardReset << 7) | (ChannelNum << 4) | CardNum;
  IOwrite(eeprom_read_word(&IOBaseAddress) - 1, ISAData);
}

void resetCardPointer(void) {
  // Resets the circular buffer data pointer for the selected card.
  if (UART_Receivei()) {
    CardReset = 1;
  } else {
    CardReset = 0;
  }
  writeCardControlRegister();
  if (verbose) {
    UART_Transmits("Card Number : ");
    UART_Transmiti(CardNum);
    UART_Transmits(", Channel Number : ");
    UART_Transmiti(ChannelNum);
    UART_Transmits(" is selected & buffer pointer is ");
    if (CardReset) {
      UART_Transmits("held at 0 ");
    } else {
      UART_Transmits("free to increment ");
    }
  }
}

void select_channel(void) {
  // Sets the selected channel.
  ChannelNum = UART_Receivei();
  writeCardControlRegister();
  if (verbose) {
    UART_Transmits("Card Number : ");
    UART_Transmiti(CardNum);
    UART_Transmits(", Channel Number : ");
    UART_Transmiti(ChannelNum);
    UART_Transmits(" is selected & buffer pointer is ");
    if (CardReset) {
      UART_Transmits("held at 0 ");
    } else {
      UART_Transmits("free to increment ");
    }
  }
}

void select_card(void) {
  // Sets the selected card.
  CardNum = UART_Receivei();
  writeCardControlRegister();
  if (verbose) {
    UART_Transmits("Card Number : ");
    UART_Transmiti(CardNum);
    UART_Transmits(", Channel Number : ");
    UART_Transmiti(ChannelNum);
    UART_Transmits(" is selected & buffer pointer is ");
    if (CardReset) {
      UART_Transmits("held at 0 ");
    } else {
      UART_Transmits("free to increment ");
    }
  }
}

void send_encoded_word(unsigned int data) {
  // The 16-bit binary data really only contains 12 bits of useful data.
  // Break this into two 6-bit values and add 0x30 to each.
  // Each character will then be contained within 7-bits and
  // be in the range of printable characters. 
  UART_Transmit((data & 0b111111) + 0x30);
  UART_Transmit(((data >> 6) & 0b111111) + 0x30);
}

unsigned int rotate_right( unsigned int c ) {
  // rotate-right with carry for 16-bit quantity
  if ( c & 1 ) 
    return (c>>1) | 0x8000;
  else
    return c>>1;
}

void getData(void) {
  // Sends the data for the selected channel on the selected card.
  // The full packet consists of a 20-byte header, a 16-kbyte payload
  // and an 8-byte trailer.
  //
  //
  if ( card_is_sampling(CardNum) ) {
    UART_Transmits("FAILED");
    return;
  }
  //
  writeCardControlRegister();
  if (verbose) {
    UART_Transmits("Card Number : ");
    UART_Transmiti(CardNum);
    UART_Transmits(", Channel Number : ");
    UART_Transmiti(ChannelNum);
    UART_Transmits(" is selected & buffer pointer is ");
    if (CardReset) {
      UART_Transmits("held at 0 ");
    } else {
      UART_Transmits("free to increment ");
    }
  }
  //
  // Header:
  // 1 Card Number = (1 -> 5)
  UART_Transmiti(CardNum);
  // 2 Channel Number = (1 -> 3)
  UART_Transmiti(ChannelNum);
  // 3 Timebase = (0 -> 3)
  unsigned char n_timebase = whichTimeBase();
  // 4 Period (us) High Byte (BCD)
  // 5 Period (us) Low Byte (BCD)
  samplePeriod(n_timebase);
  // 6 Multiplier = (0 -> 1)
  timeBaseMultiplier(n_timebase);
  // 7 Pretrigger Setting High Byte (BCD) = (0 -> 9)
  // 8 Pretrigger Setting Byte (BCD) = (0 -> 7)
  // 9 Pretrigger Setting Byte (BCD) = 9
  // 10 Pretrigger Setting Low Byte (BCD) = 2
  preTriggerDelay(n_timebase);
  // 11 Buffer Switch = (0 ,2 ,4, 8)
  bufferSizeSelect(n_timebase);
  // 12 Trigger = (1 -> 3)
  unsigned char n_trigger = whichTriggerUnit(n_timebase);
  // 13 Raw Slope and Coupling ( S C ) (0 -> 3)
  reportSlopeAndCouple(n_trigger);
  // 14 Trigger Level % ( + / - )
  // 15 Trigger Level % (tens, BCD)
  // 16 Trigger Level % (ones, BCD)
  report_threshold_in_BCD(n_trigger);
  // 17 DATA COUPLING A/D
  dataCoupling();
  // 18 FULL SCALE DATA RANGE 0 - 5
  // 19 FULL SCALE DATA RANGE .
  // 20 FULL SCALE DATA RANGE 5 - 0
  fullScaleRange();
  //
  // Payload of 16-kbytes
  // from oldest data
  //    21 B B b b b b + 0x30
  //    22 B B B B B B + 0x30
  //  ...
  // to newest data
  // 16403 B B b b b b + 0x30
  // 16404 B B B B B B + 0x30
  int buf = getBufferPointer();
  unsigned int data;
  unsigned char low_byte, high_byte;
  unsigned int i;
  unsigned int chksum = 0;
  // Start with the oldest data.
  for ( i = buf; i < 8192; i++ ) {
    low_byte = MEMread(2 * i);
    high_byte = MEMread(2 * i + 1);
    // Assemble a 12-bit integer value from the two bytes.
    data = ((unsigned int)low_byte >> 4) | ((unsigned int)high_byte << 4);
    send_encoded_word(data);
    // accumulate 16-bit checksum
    chksum = rotate_right( chksum );
    chksum += data;
    chksum &= 0xffff; // probably redundant for AVR
  }
  // Wrap around to the newest data.
  for ( i = 0; i < buf; i++ ) {
    low_byte = MEMread(2 * i);
    high_byte = MEMread(2 * i + 1);
    data = ((unsigned int)low_byte >> 4) | ((unsigned int)high_byte << 4);
    send_encoded_word(data);
    // accumulate 16-bit checksum
    chksum = rotate_right( chksum );
    chksum += data;
    chksum &= 0xffff; // probably redundant for AVR
  }
  //
  // Trailer
  UART_Transmitlong(chksum); // 16405-16408
  UART_Transmits("zzzz");    // 16409-16412 reserved
} // end getData()

unsigned int newestByte(void) {
  // Returns the most-recent word from the circular buffer.
  int buf = getBufferPointer();
  int latest;
  unsigned char high_byte, low_byte;
  unsigned int Data;
  switch (buf) {
    case 0 : latest = 16383; break;
    default : latest = (2 * buf) - 1;
  }
  high_byte = MEMread(latest);
  Data = (unsigned int)high_byte << 8;
  low_byte = MEMread(latest - 1);
  Data = (unsigned int)low_byte | Data;
  return Data;
}

void octalOutput(void) {
  // Sends an octal representation of the most-recent voltage
  // of the currently selected card and channel.
  //
  // If a card is not present (or functioning), we will most likely
  // get all high bits and an octal value of 77777.
  select_channel();
  unsigned int latest = newestByte();
  unsigned char i;
  for ( i = 5; i > 0; i--) {
    UART_Transmiti(((latest >> (3 * i - 2)) & 0b111));
  }
}

void voltageOutput(void) {
  // Sends a decimal representation of the most-recent voltage
  // of the currently selected card and channel.
  select_channel();
  unsigned int latest = newestByte();
  unsigned int mantissa = latest >> 4;
  if (mantissa > 2047) {
    UART_Transmit('+');
    mantissa = mantissa - 2048;
  } else {
    UART_Transmit('-');
    mantissa = 2048 - mantissa;
  }
  latest = (latest >> 1) & 0x7;
  mantissa = (mantissa * 31 + mantissa / 4) >> 5;
  switch (latest) {
    case 0 : mantissa = 60000; break;
    case 1 : mantissa = mantissa * 25; break;
    case 2 : mantissa = mantissa * 10; break;
    case 3 : mantissa = mantissa * 5; break;
    case 4 : mantissa = (mantissa * 5) / 2; break;
    case 5 : mantissa = mantissa * 1; break;
    case 6 :
    case 7 : mantissa = mantissa / 2; break;
  }
  UART_TransmitBCD(mantissa);
}

void fullScaleRange(void) {
  // Sends the full-scale voltage range as a fixed-point decimal number.
  unsigned int latest = (newestByte() >> 1) & 0b111 ;
  int c1[8] = {0, 5, 2, 1, 0, 0, 0, 0};
  int c2[8] = {0, 0, 0, 0, 5, 2, 1, 0};
  UART_Transmiti(c1[latest]);
  UART_Transmit('.');
  UART_Transmiti(c2[latest]);
}

void dataCoupling(void) {
  // Sends a text representation of the voltage coupling (AC/DC).
  unsigned int latest = (newestByte () & 0b1);
  if (latest) {
    UART_Transmit('D');
  } else {
    UART_Transmit('A');
  }
  if (verbose) {
    UART_Transmits("C Data Coupling");
  }
}


// ------------------- A/D Card Status ---------------------

unsigned char readCardStatusRegisterHigh(unsigned char n) {
  return IOread(eeprom_read_word(&IOBaseAddress) + 2 * n + 1);
}

unsigned char readCardStatusRegisterLow(unsigned char n) {
  return IOread(eeprom_read_word(&IOBaseAddress) + 2 * n);
}

unsigned char card_is_present(unsigned char n) {
  // Returns 1 if card n appears to be present, 0 otherwise.
  // If a card is missing, the ISA bus terminators will pull
  // the data lines high.
  unsigned char all_ones = (readCardStatusRegisterHigh(n) == 0xFF) &&
                           (readCardStatusRegisterLow(n) == 0xFF);
  return !all_ones;
}

unsigned char card_is_sampling(unsigned char n) {
  // Returns 1 is still sampling, 0 otherwise.
  unsigned char ISAData = readCardStatusRegisterHigh(n);
  return ISAData >> 7;
}

unsigned char any_cards_sampling(void) {
  // Returns 1 if any existing cards are sampling, 0 otherwise.
  unsigned char any_sampling = 0;
  unsigned char i;
  for ( i = 1; i <= 7; ++i ) {
    if ( card_is_present(i) ) {
      if ( card_is_sampling(i) ) any_sampling = 1;
    }
  }
  return any_sampling;
}

void report_sampling_status(unsigned char n) {
  // Sends text indicating if the card is still sampling.
  if (verbose) {
    UART_Transmits("Card Number : ");
    UART_Transmiti(n);
    UART_Transmits(" sampling : ");
  }
  UART_Transmiti(card_is_sampling(n));
}

unsigned char whichTimeBase(void) {
  // Works out which timebase is being used by the currently-selected card.
  // Sends a text representation of the timebase number.
  unsigned char ISAData = readCardStatusRegisterHigh(CardNum);
  int n_timebase = (ISAData >> 5) & 0b11;
  if (verbose) {
    UART_Transmits("Card Number : ");
    UART_Transmiti(CardNum);
    UART_Transmits(" is using timebase ");
  }
  UART_Transmiti(n_timebase);
  return n_timebase;
}

int getBufferPointer (void) {
  // Gets the value of the poiner to the oldest word in the circular buffer.
  int buf;
  unsigned char ISAData = readCardStatusRegisterHigh(CardNum);
  buf = (ISAData & 0b11111) << 8;
  ISAData = readCardStatusRegisterLow(CardNum);
  buf = buf | ISAData;
  return buf;
}

void printBufferPointer(void) {
  // Returns a text representation of the pointer to the oldest word
  // in the circular buffer.
  if (verbose) {
    UART_Transmits("Card Number : ");
    UART_Transmiti(CardNum);
    UART_Transmits(" buffer pointer is 0x");
  }
  UART_Transmitlong(getBufferPointer());
}


// ---------------- Trigger Unit Control ------------------

void writeTriggerUnitControlRegister(unsigned char ISAData) {
  IOwrite(eeprom_read_word(&TimeBaseAddress), ISAData);
  return;
}

void arm_box(void) {
  // Start the box sampling by writing the arm bit to the register.
  if (UART_Receivei()) {
    writeTriggerUnitControlRegister( 0x20 );
    if (verbose) UART_Transmits("All Trigger Units Armed.");
  }
}

void trigger_box(void) {
  // Stop the box sampling (sometime soon) by writing the trigger bit.
  if (UART_Receivei()) {
    writeTriggerUnitControlRegister( 0x10 );
    if (verbose) UART_Transmits("All Trigger Units Triggered.");
  }
}

// ------------------- Time Base Status ------------------------

unsigned char readTimeStatusRegisterHigh(unsigned char n) {
  if (n <= 0) n = 1;
  if (n > 3) n = 3;
  return IOread(eeprom_read_word(&TimeBaseAddress) + 1 + 2 * n);
}

unsigned char readTimeStatusRegisterLow(unsigned char n) {
  if (n <= 0) n = 1;
  if (n > 3) n = 3;
  return IOread(eeprom_read_word(&TimeBaseAddress) + 2 * n);
}

unsigned char box_is_present() {
  // Returns 1 if box appears to be present, 0 otherwise.
  unsigned char all_ones = (readTimeStatusRegisterHigh(1) == 0xFF) &&
                           (readTimeStatusRegisterLow(1)  == 0xFF) &&
                           (readTimeStatusRegisterHigh(2) == 0xFF) &&
                           (readTimeStatusRegisterLow(2)  == 0xFF) &&
                           (readTimeStatusRegisterHigh(3) == 0xFF) &&
                           (readTimeStatusRegisterLow(3)  == 0xFF);
  return !all_ones;
}

void preTriggerDelay(unsigned char n) {
  // Sends a text representation of the number of pretrigger samples
  // for timebase n.
  if (n <= 0) n = 1;
  if (n > 3) n = 3;
  unsigned char ISAData = readTimeStatusRegisterHigh(n);
  if (verbose) {
    UART_Transmits("Trigger Unit ");
    UART_Transmiti(n);
    UART_Transmits(" Pre-trigger samples is ");
  }
  UART_Transmiti((ISAData >> 4) & 0b111);
  UART_Transmiti(ISAData & 0b1111);
  UART_Transmits("92");
}

void bufferSizeSelect(unsigned char n) {
  // Sends a text representation of the selected buffer size for timebase n.
  if (n <= 0) n = 1;
  if (n > 3) n = 3;
  unsigned char ISAData = readTimeStatusRegisterLow(n);
  if (verbose) {
    UART_Transmits("Time Base Unit ");
    UART_Transmiti(n);
    UART_Transmits(" buffer size is ");
  }
  switch (ISAData >> 6) {
    case 0 : UART_Transmits("A"); break;
    case 1 : UART_Transmits("2"); break;
    case 2 : UART_Transmits("8"); break;
    case 3 : UART_Transmits("4"); break;
    default : UART_Transmits("X"); break;
  } // end switch
}

unsigned char whichTriggerUnit(unsigned char n) {
  // Send a text representation of the trigger-unit that is being used
  // by the timebase n.
  unsigned char trigger_unit;
  if (n <= 0) n = 1;
  if (n > 3) n = 3;
  unsigned char ISAData = readTimeStatusRegisterLow(n);
  if (verbose) {
    UART_Transmits("Time Base Unit ");
    UART_Transmiti(n);
    UART_Transmits(" is using trigger unit ");
  }
  switch ((ISAData >> 4) & 0b11) {
    // Note the odd mapping; it's not straight through.
    case 0 : trigger_unit = 1; break;
    case 1 : trigger_unit = 2; break;
    case 3 : trigger_unit = 3; break;
    default : trigger_unit = 0; break; // should never see this case
  }
  UART_Transmiti(trigger_unit);
  return trigger_unit;
}

void samplePeriod(unsigned char n) {
  // Sends a text representation of the sample period (in us)
  // for the timebase n.
  if (n <= 0) n = 1;
  if (n > 3) n = 3;
  unsigned char ISAData = readTimeStatusRegisterLow(n);
  if (verbose) {
    UART_Transmits("Trigger Unit ");
    UART_Transmiti(n);
    UART_Transmits(" sampling period is (us) ");
  }
  switch (ISAData & 0b111) {
    case 1 : UART_Transmits("50"); break;
    case 2 : UART_Transmits("20"); break;
    case 3 : UART_Transmits("10"); break;
    case 4 : UART_Transmits("05"); break;
    case 5 : UART_Transmits("02"); break;
    case 6 : UART_Transmits("01"); break;
    default : UART_Transmits("XX"); break;
  }
}

void timeBaseMultiplier(unsigned char n) {
  // Sends a text representation of the sample-period multiplier
  // for timebase n.
  if (n <= 0) n = 1;
  if (n > 3) n = 3;
  unsigned char ISAData = readTimeStatusRegisterLow(n);
  if (verbose) {
    UART_Transmits("Trigger Unit ");
    UART_Transmiti(n);
    UART_Transmits(" time base multiplier is  ");
  }
  switch ((ISAData >> 3) & 0b1) {
    case 0 : 
      UART_Transmits("0"); break;  //x1 multiplier
    case 1 : 
      UART_Transmits("1");       //x100 multiplier
      if (verbose) {
        UART_Transmits(" x100");
      }
      break;
    default : UART_Transmits("X"); break;
  }
}

// ------------------ Trigger Unit Status Registers --------------------------

unsigned char readTriggerStatusRegisterHigh(void) {
  return IOread(eeprom_read_word(&TimeBaseAddress) + 1);
}

unsigned char readTriggerStatusRegisterLow(void) {
  return IOread(eeprom_read_word(&TimeBaseAddress));
}

void report_coupling(unsigned char n) {
  // Sends a text representation of the coupling for trigger unit n.
  if ( n < 1 ) n = 1;
  if ( n > 3 ) n = 3;
  unsigned char ISAData = readTriggerStatusRegisterLow();
  if (verbose) {
    UART_Transmits("Trigger : ");
    UART_Transmiti(n);
    UART_Transmit(' ');
  }
  if ((ISAData >> (2 * n - 2)) & 0b1) {
    UART_Transmit('D');
  } else {
    UART_Transmit('A');
  }
  if (verbose) {
    UART_Transmits("C Coupling");
  }
}

void report_slope(unsigned char n) {
  // Sends a text representation of the slope for trigger unit n.
  if ( n < 1 ) n = 1;
  if ( n > 3 ) n = 3;
  unsigned char ISAData = readTriggerStatusRegisterLow();
  if (verbose) {
    UART_Transmits("Slope on Trigger : ");
    UART_Transmiti(n);
    UART_Transmits(" ");
  }
  if((ISAData >> (2 * n - 1)) & 0b1) {
    UART_Transmits("R"); // rising
  } else {
    UART_Transmits("F"); // falling
  }
}

void reportSlopeAndCouple(unsigned char n) {
  // Sends a text representation of the raw bits for slope and coupling
  // of trigger unit n.
  if ( n < 1 ) n = 1;
  if ( n > 3 ) n = 3;
  unsigned char ISAData = readTriggerStatusRegisterLow();
  if (verbose) {
    UART_Transmits("Slope and Coupling on Trigger : ");
    UART_Transmiti(n);
    UART_Transmits(" ");
  }
  UART_Transmiti((ISAData >> (2 * n - 2)) & 0b11);
}

void enable_threshold_report(unsigned char n) {
  // Set the flag for suppressing or allowing reporting of k.
  // I suppose that the logic looks a bit backward but
  // it is easier to short-cut the reporting functions if
  // k-reporting is to be suppressed.
  if ( n == 1 ) {
    suppress_k_report = 0;
  } else {
    suppress_k_report = 1;
  }
  if (verbose) {
    UART_Transmits("Enable k-reporting: ");
    UART_Transmiti(!suppress_k_report);
    UART_Transmits(" ");
  } else {
    UART_Transmiti(!suppress_k_report);
  }
  return;
}

void report_threshold(unsigned char n) {
  // Sends an octal representation of the threshold for trigger unit n.
  //
  // There is a small catch with the databox: it will trigger if we
  // write to the control register to select the requested unit.
  if ( any_cards_sampling() || suppress_k_report ) {
    UART_Transmits("000");  // a null value, so to speak
    return;
  }
  if ( n < 1 ) n = 1;
  if ( n > 3 ) n = 3;
  unsigned char ISAData = (0x04 | n) & 0x0F; // redundant filtering of bits
  writeTriggerUnitControlRegister(ISAData);
  while ( thresholdStatus() ) /* wait */ ;
  ISAData = readTriggerStatusRegisterHigh();
  if (verbose) {
    UART_Transmits("Threshold on Trigger : ");
    UART_Transmiti(n);
    UART_Transmits(" ");
  }
  // OCTAL 000 = -2.5V, 377 = 2.5V
  char i;
  for (i = 2; i >= 0; i--) {
    UART_Transmit(((ISAData >> (3 * i)) & 0b111) + 0x30);
  }
}

void report_threshold_in_BCD(unsigned char n) {
  // Sends a decimal representation of the threshold for trigger unit n.
  //
  // There is a small catch with the databox: it will trigger if we
  // write to the control register to select the requested unit.
  if ( any_cards_sampling() || suppress_k_report ) {
    UART_Transmits("+00");  // a null value, so to speak
    return;
  }
  if ( n < 1 ) n = 1;
  if ( n > 3 ) n = 3;
  unsigned char ISAData = (0x04 | n) & 0x0F;
  writeTriggerUnitControlRegister(ISAData);
  while ( thresholdStatus() ) /* wait */ ;
  ISAData = readTriggerStatusRegisterHigh();
  if (verbose) {
    UART_Transmits("Threshold on Trigger : ");
    UART_Transmiti(n);
    UART_Transmits(" ");
  }
  //-99 = -2.5V, +99 = 2.5V
  if (ISAData > 127) {
    UART_Transmit('+');
    ISAData = ISAData - 128;
  } else {
    UART_Transmit('-');
    ISAData = 127 - ISAData;
  }
  unsigned int value = ISAData * 100 / 128;
  unsigned int high = value / 10;
  UART_Transmiti(high);
  UART_Transmiti(value - high * 10);
}

unsigned char thresholdStatus(void) {
  // Returns 1 if the threshold A->D converter is busy.
  unsigned char ISAData = readTriggerStatusRegisterLow();
  return (ISAData >> 7);
}

// ----------------------- System Functions ------------------------

void reportIOBase(void) {
  // Sends a text representation of the current IO base-address.
  if (verbose) {
    UART_Transmits("I/O BASE ADDRESS : 0x");
  }
  UART_Transmitlong(eeprom_read_word(&IOBaseAddress));
}

void setIOBase(void) {
  // Sets the base-addresses in the IO space by writing them into the EEPROM.
  int temp;
  temp =  UART_Receivei() << 8;
  temp |= UART_Receivei() << 4;
  temp |= UART_Receivei();
  eeprom_write_word(&IOBaseAddress, temp);
  eeprom_write_word(&TimeBaseAddress, (temp - 0x10));
  reportIOBase();
}

void setVerbose(void) {
  // Sets the status of the verbose flag.
  verbose = UART_Receivei();
  if (verbose) {
    UART_Transmits("Verbose : ON");
  } else {
    UART_Transmit('F');
  }
}

void menu(void) {
  const char *menu_text =
#   include "Menu.txt"
  ;
  char c;
  const char *addr;
  addr = menu_text;
  while ( (c = pgm_read_byte(addr++)) ) {
    UART_Transmit(c);
  }
}
