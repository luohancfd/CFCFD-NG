/* UART-INTERRUPT.C
 * An implmentation of a interrupt driven UART for the Databox Interface
 * Project. Uses circular buffers for RX and TX and automatically merges data
 * from both UARTS.
 *
 * 15/10/2006 Luke Hillyard
 * 16-Jan-2007 Peter J.  clean-up and comment
 */

#include <avr/interrupt.h>
#include <avr/io.h>
#include <avr/eeprom.h>
#include "bits.h"
#include "databox.h"
#include "uart-interrupt.h"

//Constants for Baud Rate Calculation
extern unsigned char source;

//Size of RX and TX buffers in bytes (chars)
#define Uart0TXBuffSizeMax 64
#define Uart0RXBuffSizeMax 64
#define Uart1TXBuffSizeMax 64
#define Uart1RXBuffSizeMax 64

//RX and TX buffers;
volatile char Uart0TXBuff[Uart0TXBuffSizeMax];
volatile char Uart0RXBuff[Uart0RXBuffSizeMax];
volatile char Uart1TXBuff[Uart1TXBuffSizeMax];
volatile char Uart1RXBuff[Uart1RXBuffSizeMax];

//Variables to track the size of the buffers
volatile unsigned char Uart0TXBuffSize = 0;
volatile unsigned char Uart0RXBuffSize = 0;
volatile unsigned char Uart1TXBuffSize = 0;
volatile unsigned char Uart1RXBuffSize = 0;

//Pointer to Unread / Unsend data
volatile unsigned char Uart0TXBuffLoc = 0;
volatile unsigned char Uart0RXBuffLoc = 0;
volatile unsigned char Uart1TXBuffLoc = 0;
volatile unsigned char Uart1RXBuffLoc = 0;

//Pointer to next free location in buffer
volatile unsigned char Uart0TXBuffInsertLoc = 0;
volatile unsigned char Uart0RXBuffInsertLoc = 0;
volatile unsigned char Uart1TXBuffInsertLoc = 0;
volatile unsigned char Uart1RXBuffInsertLoc = 0;


// -------------- Interrupt Handlers -------------------

SIGNAL (SIG_UART0_RECV) {
  // A character has been received from the USB port; put it in the receive buffer.
  cli();
  if (Uart0RXBuffSize < Uart0RXBuffSizeMax) { 
    // Space is available
    Uart0RXBuffSize++;
    Uart0RXBuff[Uart0RXBuffInsertLoc] = UDR0;
    Uart0RXBuffInsertLoc++;
    if (Uart0RXBuffInsertLoc == Uart0RXBuffSizeMax) {
      Uart0RXBuffInsertLoc = 0;
    }
  } else {
    /* do nothing */ ; // Lose the character
  }
  sei();
}

SIGNAL (SIG_UART0_DATA) {
  // The transmit side of UART0 is waiting for more data.
  // Send a character if there is one in the buffer, 
  // else tell UART0 to stop interrupting.
  cli();
  // Presumably, we have data to send.
  UDR0 = Uart0TXBuff[Uart0TXBuffLoc];
  Uart0TXBuffSize--;
  Uart0TXBuffLoc++;
  if (Uart0TXBuffLoc == Uart0TXBuffSizeMax) {
    Uart0TXBuffLoc = 0;
  }
  if (Uart0TXBuffSize == 0) { 
    // There is no more data to send; tell UART0 to stop interrupting.
    CLEARBIT(UCSR0B,UDRIE);
  }
  sei();
}

SIGNAL (SIG_UART1_RECV) {
  // A character has been received at UART1; put it in the receive buffer.
  cli();
  if (Uart1RXBuffSize < Uart1RXBuffSizeMax) { 
    // There is space available.
    Uart1RXBuff[Uart1RXBuffInsertLoc] = UDR1;
    Uart1RXBuffSize++;
    Uart1RXBuffInsertLoc++;
    if (Uart1RXBuffInsertLoc == Uart1RXBuffSizeMax) {
      Uart1RXBuffInsertLoc = 0;
    }
  } else {
    /* do nothing */ ; // Lose the character.
  }
  sei();
}

SIGNAL (SIG_UART1_DATA) {
  // The transmit side of UART1 is waiting for more data.
  // Send a character if there is one in the buffer, 
  // else tell UART1 to stop interrupting.
  cli();
  // Presumably, we have data to send.
  UDR1 = Uart1TXBuff[Uart1TXBuffLoc];
  Uart1TXBuffSize--;
  Uart1TXBuffLoc++;
  if (Uart1TXBuffLoc == Uart1TXBuffSizeMax) {
    Uart1TXBuffLoc = 0;
  }
  if (Uart1TXBuffSize == 0) { 
    // No data to send; tell UART1 to stop interrupting.
    CLEARBIT(UCSR1B,UDRIE);
  }
  sei();
}

// ------------- functions that deal directly with hardware -------------

void UART_Init (unsigned char rs232_baud_rate_selection) {
  // The USART Initialisation routine
  //
  // Both UARTs:
  //   Asynchronous Normal Mode
  //   7 Data Bits, Odd Parity, 2 Stop Bits
  //   No Flow Control
  // Set the baud rates for each of the UARTs:
  #define FOSC 18432000L
  unsigned int UBRR;
  unsigned long rs232_baud_rate;

  cli(); // we don't want to be interrupted while fiddling UART config bits.
  // UART0 (connected to USB) at 230.4k.
  #define BAUD 230400
  UBRR0H = (unsigned char)((FOSC/16/BAUD-1)>>8);
  UBRR0L = (unsigned char)(FOSC/16/BAUD-1);
  /* 7 Bit mode */
  CLEARBIT(UCSR0B, UCSZ02);
  SETBIT(UCSR0C, UCSZ01);
  CLEARBIT(UCSR0C, UCSZ00);
  /* 2 Stop bit */
  SETBIT(UCSR0C, USBS0);
  /* Odd Parity */
  SETBIT(UCSR0C, UPM01);
  SETBIT(UCSR0C, UPM00);
  /* Enable receiver and transmitter */
  SETBIT(UCSR0B,RXEN);  //Enable Recieve
  SETBIT(UCSR0B,TXEN);  //Enable Transmit
  SETBIT(UCSR0B,RXCIE); //Interupt when data recieved

  // UART1 (connected via RS232) baud-rate encoded in UBRR.
  switch (rs232_baud_rate_selection) {
    case 1: rs232_baud_rate = 2400; break;
    case 2: rs232_baud_rate = 9600; break;
    case 3: rs232_baud_rate = 19200; break;
    case 4: rs232_baud_rate = 38400; break;
    case 5: rs232_baud_rate = 57600; break;
    case 6: rs232_baud_rate = 115200L; break;
    default: rs232_baud_rate = 57600;
  }
  UBRR = (FOSC/16/rs232_baud_rate - 1);

  UBRR1H = (unsigned char)(UBRR>>8);
  UBRR1L = (unsigned char)UBRR;
  /* 7 Bit mode */
  CLEARBIT(UCSR1B, UCSZ12);
  SETBIT(UCSR1C, UCSZ11);
  CLEARBIT(UCSR1C, UCSZ10);
  /* 2 Stop bit */
  SETBIT(UCSR1C, USBS1);
  /* Odd Parity */
  SETBIT(UCSR1C, UPM11);
  SETBIT(UCSR1C, UPM10);
  /* Enable receiver and transmitter */
  SETBIT(UCSR1B,RXEN);  //Enable Recieve
  SETBIT(UCSR1B,TXEN);  //Enable Transmit
  SETBIT(UCSR1B,RXCIE); //Interupt when data recieved

  sei(); // let the interrupt routines take over.
} // end UART_Init()


void UART_Transmit0 (unsigned char data) {
  // Send a character to the USB port circular buffer.
  while (Uart0TXBuffSize >= Uart0TXBuffSizeMax) { 
    // The TX buffer is full; let the hardware deal with it.
    ; // Wait
  }
  cli();
  Uart0TXBuff[Uart0TXBuffInsertLoc] = data;
  Uart0TXBuffSize++;
  Uart0TXBuffInsertLoc++;
  if (Uart0TXBuffInsertLoc == Uart0TXBuffSizeMax) {
    Uart0TXBuffInsertLoc = 0;
  }
  // Since there is now data to send, tell UART1 to signal when it is ready.
  SETBIT(UCSR0B,UDRIE);
  sei();
}

void UART_Transmit1 (unsigned char data) {
  // Transmit a character to the RS232 serial port.
  while (Uart1TXBuffSize >= Uart1TXBuffSizeMax) { 
    // The TX buffer is full; let the hardware deal with it.
    ; // Wait
  }
  cli();
  Uart1TXBuff[Uart1TXBuffInsertLoc] = data;
  Uart1TXBuffSize++;
  Uart1TXBuffInsertLoc++;
  if (Uart1TXBuffInsertLoc == Uart1TXBuffSizeMax) {
    Uart1TXBuffInsertLoc = 0;
  }
  // Since there is now data to send, tell UART1 to signal when it is ready.
  SETBIT(UCSR1B,UDRIE);
  sei();
}

unsigned char UART_Receive0 (void) {
  // Receive a character from the USB port.
  unsigned char data;
  while (Uart0RXBuffSize == 0) { 
    // The receive buffer is empty; wait for the hardware to put something in it.
    /* do nothing */ ;
  }
  cli();
  data = Uart0RXBuff[Uart0RXBuffLoc];
  Uart0RXBuffSize--;
  Uart0RXBuffLoc++;
  if (Uart0RXBuffLoc == Uart0RXBuffSizeMax) {
    Uart0RXBuffLoc = 0;
  }
  sei();
  return data;
}

unsigned char UART_Receive1 (void) {
  // Receive a character from the RS232 port.
  unsigned char data;
  while (!Uart1RXBuffSize) {
    // The receive buffer is empty; wait for the hardware to put something in it.
    /* do nothing */ ;
  }
  cli();
  data = Uart1RXBuff[Uart1RXBuffLoc];
  Uart1RXBuffSize--;
  Uart1RXBuffLoc++;
  if (Uart1RXBuffLoc == Uart1RXBuffSizeMax) {
    Uart1RXBuffLoc = 0;
  }
  sei();
  return data;
}

// ----------- functions to unify the UART streams --------------

void UART_Transmit (unsigned char data) {
  // Send a character to the appointed serial port.
  if (source) {
    UART_Transmit1(data);
  } else {
    UART_Transmit0(data);
  }
}

unsigned char UART_Receive (void) {
  // Hang around until we get one character from either UART.
  for (;;) {
    if (Uart0RXBuffSize) {
      source = 0;
      return UART_Receive0();
    }
    if (Uart1RXBuffSize) {
      source = 1;
      return UART_Receive1();
    }
  } // end for
}

// ------ higher-level functions to deal with various data types ------

void UART_TransmitBCD (unsigned int V) {
  // Send a 16-bit fixed-point value as a decimal digits, 
  // starting with the most significant.
  unsigned int i;
  unsigned char C;
  for (i = 10000; i >= 1; i = i / 10) {
    C = V / i;
    V = V - C * i;
    UART_Transmiti(C);
    if (i == 10000) {
      UART_Transmit('.');
    }
  }
}

void UART_Transmiti (unsigned char data) {
  // Send a binary value encoded as a single hexadecimal digit.
  unsigned char temp;
  temp = (unsigned char)(data & 0xF);
  if (temp < 10) {
    temp = temp + '0';
  } else {
    temp = temp + 'A' - 10;
  }
  UART_Transmit(temp);
}

void UART_Transmitlong (unsigned int data) {
  // Send a 16-bit integer as 4 hexadecimal digits, one for each nibble,
  // starting with the most significant.
  char temp;
  char i;
  for (i = 3; i >= 0; i--) {
    temp = (data >> (4*i)) & 0xF;
    if (temp < 10) {
      temp = temp + '0';
    } else {
      temp = temp + 'A' - 10;
    }
    UART_Transmit(temp);
  }
}

void UART_Transmits (unsigned char* data) {
  // Send a null-terminated string (but not the null character itself).
  while (*data != '\0') {
    UART_Transmit(*data++);
  }
}

void UART_Transmitb (unsigned char data) {
  // Send a byte, one bit at a time, encoded as ASCII '1' or '0'.
  char i;
  for (i = 7; i >= 0; i--) {
    UART_Transmit(((data >> i) & 0b1) + '0');
  }
}

int UART_Receivei (void) {
  // Retrieve one hexadecimal digit from the UART and return its binary value.
  unsigned char temp = UART_Receive();
  if (temp >= 'A') {
    return (int)((temp - 'A' + 10) & 0xF);
  } else {
    return (int)((temp - '0') & 0xF);
  }
}

