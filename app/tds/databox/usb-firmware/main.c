/* main.c
 * A project to provide an interface between the ISA bus and the Serial and
 * USB ports. Operating at the highest level this file maps characters
 * received from the serial port to functions that return data as required.
 *
 * 15/10/2006 Luke Hillyard
 * 16-Jan-2007 Peter J. -- code clean-up and comment
 * 17-Jan-2007 Peter J. -- some real restructuring in all of the files
 * 19-Jan-2007 Peter J. -- finalize and add checksum calculation
 * 24-Jan-2007 Peter J. -- rework assembly of data values from memory
 * 29-Jan-2007 Peter J. -- v3.3 Suppress reporting of threshold for T4 databox.
 * 05-Feb-2007 Peter J. -- v3.4 command to enable/disable threshold-reporting
 * 06-Feb-2007 Peter J. -- v3.5 move the menu text to program space so that we don't
 *                         consume all of the RAM with string data.
 */

#define VERSION "3.5"

#include <avr/io.h>
#include <avr/eeprom.h>
#include "isa-bus.h"
#include "databox.h"
#include "uart-interrupt.h"
#include "main.h"


unsigned char EEMEM baudRateSelection = 6; // 115200 baud
unsigned char verbose = 1; // Verbose = 1
unsigned char source = 1;  // RS232 Serial = 1
                           // USB Serial = 0

void setBaudRate () {
  // Reset the RS232 baud-rate from the selection received.
  eeprom_write_byte(&baudRateSelection, UART_Receivei());
  UART_Init(eeprom_read_byte(&baudRateSelection));
}

// --------------------------------------------------------------------------------

int main (void) {
  ISA_Init();
  UART_Init(eeprom_read_byte(&baudRateSelection));

  // The program is now driven by the data that comes down the serial line.
  for (;;) {
    switch ( UART_Receive() ) {
      // All uppercase and some lowercase commands require aninteger parameter (i)
      // e.g. R1 to reset //

      // A/D Card Control
      case 'R' : resetCardPointer(); break;  //<i=1-0> 0=reset & hold ring pointers at zero.
      case 'N' : select_card(); break;  //<i=1-F> select card for other commands
      case 'C' : select_channel(); break;  //<i=1-3> select channel for other commands

      case 'D' : verbose = 0;   //force brief mode for all subsequent commands
             select_channel();
             getData(); break;  // send header+8k values for card N

      case 'O' : octalOutput(); break;  //<i=1-3> transmit the latest 2 bytes as octal data
      case 'V' : voltageOutput(); break;  //<i=1-3> transmit the latest 2 bytes as BCD data
      case 'f' : fullScaleRange(); break; //transmit the latest data range
      case 'g' : dataCoupling(); break;  //transmit the latest data coupling

      // A/D Card Status
      case 'x' : UART_Transmiti(card_is_present(UART_Receivei())); break; // sends '1' if card is present, '0' otherwise
      case 'a' : report_sampling_status(UART_Receivei()); break;  // 1=still sampling, 0=hold
      case 't' : whichTimeBase(); break;    //switch: 1 to 3 or 0=timebase#1 x4
      case 'r' : printBufferPointer(); break;  //get Ring buffer oldest data address,0xBBBB

      // Trigger Unit Control
      case 'A' : arm_box(); break;             //<i=1> 1=arm all trigger units
      case 'T' : trigger_box(); break;    //<i=1> 1=causes a triggering of databox

      // Time Base Status
      case 'd' : preTriggerDelay(UART_Receivei()); break;  //get thumbwheels for time base 0xCD92
      case 'b' : bufferSizeSelect(UART_Receivei()); break;  //get selected size (2,4,8kb,A=Active)
      case 'p' : samplePeriod(UART_Receivei()); break;    //get sample period for timebase 01-50
      case 'm' : timeBaseMultiplier(UART_Receivei());break;  //get multiplier 1=x100 0=x1

      case 'u' : whichTriggerUnit(UART_Receivei()); break;

      // Trigger Unit Status
      case 'c' : report_coupling(UART_Receivei()); break;    // D=1=dc,A=0=ac
      case 's' : report_slope(UART_Receivei()); break;       // 1=rising,0=falling
      case 'L' : report_threshold(UART_Receivei()); break;   // in octal
      case 'k' : report_threshold_in_BCD(UART_Receivei()); break; // decimal percentage  +/- 99
      case 'z' : enable_threshold_report(UART_Receivei()); break;    // 1=enable,0=disable

      // System Functions
      case 'y' : UART_Transmiti(box_is_present()); break; // sends '1' if databox is present, '0' otherwise
      case 'v' :
        if ( source && verbose ) {
          UART_Transmits("Data Box V" VERSION " Present on Serial Port");
        } else if ( source ) {
          UART_Transmits("S" VERSION);
        } else if ( verbose ) {
          UART_Transmits("Data Box V" VERSION " Present on USB Port");
        } else {
          UART_Transmits("U" VERSION);
        }
        if ( verbose ) {
          UART_Transmits("\r\n");
          UART_Transmits("Luke Hillyard: hardware and firmware, 2006.\r\n");
          UART_Transmits("Peter Jacobs: firmware, Jan 2007.\r\n");
        }
        break;
      case 'i' : reportIOBase(); break; //get programmed i/o base address: default=0x0320
      case 'I' : setIOBase(); break;    //set i/o base address. Needs 3 characters! e.g 320
                        // NOTE : BASE ADDRESS IS REMEMBERED AFTER POWER CYCLE
      case 'U' : setBaudRate(); break;  // <i=1-6> set RS232 baud rate
                        // NOTE : BAUD RATE IS REMEMBERED AFTER POWER CYCLE
      case 'B' : setVerbose(); break;   // <i=0-1> 0=brief 1=verbose
      case 'h' : menu(); break;

      default:
        if ( verbose ) {
          UART_Transmits("What the fxck");
        }
        UART_Transmit('?');
        break;
    } // end switch()

    if (verbose) {
      UART_Transmits("\r\n");
    }
  } // end for()
  return 0;
}

