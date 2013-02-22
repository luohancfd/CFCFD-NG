/* uart-interrupt.h */

#ifndef __UART_INTERRPUT_H
#define __UART_INTERRUPT_H

void UART_Init (unsigned char rs232_baud_rate_selection);
void UART_Transmit0 (unsigned char data);
void UART_Transmit1 (unsigned char data);
unsigned char UART_Receive0 (void);
unsigned char UART_Receive1 (void);

void UART_Transmit (unsigned char data);
unsigned char UART_Receive (void);

void UART_TransmitBCD (unsigned int data);
void UART_Transmiti (unsigned char data);
void UART_Transmitlong (unsigned int data);
void UART_Transmits (unsigned char* data);
void UART_Transmitb (unsigned char data);

int UART_Receivei (void);

#endif
