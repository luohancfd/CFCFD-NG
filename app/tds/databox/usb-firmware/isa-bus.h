// isa-bus.h

// Low Level Functions
void ISA_Init (void);
unsigned char IOread(unsigned int addr);
void IOwrite(unsigned int addr, unsigned char ISAData);
unsigned char MEMread(unsigned int addr);
