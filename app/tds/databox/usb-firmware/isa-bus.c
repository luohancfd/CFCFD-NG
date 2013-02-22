// isa-bus.c
// Low level functions that fiddle the bits on the ISA bus of the databox.
//
// Luke Hillyard, October 2006
// Peter Jacobs, January 2007
//

#include <avr/io.h>
#include <avr/eeprom.h>
#include "bits.h"
#include "isa-bus.h"

// Addressing memory across the ISA bus.
                        // X1 (0 = on; 1 = off)
// 0xD8000 Databox -->> 0b 1101 1000 0000 0000 0000
// Interface Card  -->> 0b xxxx xx00 0000 0000 0000
// With the way that the hardware is built, 
// there is now no way to select anything other than
// 0b 1101 10xx xxxx xxxx
// 0x D    8    0    0

#define DATABUSDIR DDRF
#define INPUT 0x00
#define OUTPUT 0xFF
#define IORC PIN1
#define IOWC PIN0
#define SMRDC PIN3

void ISA_Init(void) {
  DDRA = 0xFF;       //Address Bus Low-order bits
  // Use line below to mask out SA14 and SA15 (PINC6 and PINC7)
  DDRC = 0b00111111; //Address Bus High-order bits
  //
  // Control Lines - HIGH is the ISA default
  SETBIT(PORTG, SMRDC);
  SETBIT(DDRG, SMRDC);
  SETBIT(PORTG, IORC);
  SETBIT(DDRG, IORC);
  SETBIT(PORTG, IOWC);
  SETBIT(DDRG, IOWC);
  return;
}

unsigned char IOread(unsigned int addr) {
  // Return a byte from the ISA bus, IO area.
  //
  unsigned char ISAData;
  //
  // Set up the address pattern for the subset of the ISA bus
  // that we are using in the databox.
  PORTA = addr & 0xFF;
  PORTC = (addr>>8) & 0xFF;
  //
  DATABUSDIR = INPUT;
  CLEARBIT(PORTG,IORC);
  // Insert a 0.325us wait.
  // This should be long enough for Barry's electronics to have settled.
  asm volatile("NOP"::); // at 18.432MHz, each NOP represents 54.2ns
  asm volatile("NOP"::);
  asm volatile("NOP"::);
  asm volatile("NOP"::);
  asm volatile("NOP"::);
  asm volatile("NOP"::);
  ISAData = PINF;
  SETBIT(PORTG,IORC);
  return ISAData;
}

void IOwrite(unsigned int addr, unsigned char ISAData) {
  // Write a byte (ISAData) onto the ISA bus as IO data.
  //
  // Set up the address pattern for the subset of the ISA bus
  // that we are using in the databox.
  PORTA = addr & 0xFF;
  PORTC = (addr>>8) & 0xFF;
  //
  PORTF = ISAData;
  DATABUSDIR = OUTPUT;
  CLEARBIT(PORTG,IOWC);
  asm volatile("NOP"::);
  asm volatile("NOP"::);
  asm volatile("NOP"::);
  asm volatile("NOP"::);
  asm volatile("NOP"::);
  asm volatile("NOP"::);
  SETBIT(PORTG,IOWC);
  // PORTG = 0xFF; // Faster and the default is high anyway.
  DATABUSDIR = INPUT;
  return;
}

unsigned char MEMread(unsigned int addr) {
  // Returns a byte from the ISA memory area.
  //
  unsigned char ISAData;
  //
  // Set up the address pattern for the subset of the ISA bus
  // that we are using in the databox.
  PORTA = addr & 0xFF;
  PORTC = (addr>>8) & 0xFF;
  //
  DATABUSDIR = INPUT;
  CLEARBIT(PORTG,SMRDC);
  asm volatile("NOP"::);
  asm volatile("NOP"::);
  asm volatile("NOP"::);
  asm volatile("NOP"::);
  asm volatile("NOP"::);
  asm volatile("NOP"::);
  ISAData = PINF;
  SETBIT(PORTG,SMRDC);
  return ISAData;
}

