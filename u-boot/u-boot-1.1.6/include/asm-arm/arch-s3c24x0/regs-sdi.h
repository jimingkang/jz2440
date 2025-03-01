/* linux/include/asm/arch-s3c2410/regs-sdi.h
*
* Copyright (c) 2004 Simtec Electronics <linux@simtec.co.uk>
*                             http://www.simtec.co.uk/products/SWLINUX/ *
* This program is free software; you can redistribute it and/or modify
* it under the terms of the GNU General Public License version 2 as
* published by the Free Software Foundation. *
* S3C2410 MMC/SDIO register definitions *
*  Changelog:
*    18-Aug-2004 Ben Dooks      Created initial file
*    29-Nov-2004 Koen Martens  Added some missing defines, fixed duplicates
*    29-Nov-2004 Ben Dooks       Updated Koen's patch
*/
#ifndef __ASM_ARM_REGS_SDI
#define __ASM_ARM_REGS_SDI "regs-sdi.h"
#define S3C2440_SDICON_SDRESET        (1<<8)
#define S3C2440_SDICON_MMCCLOCK       (1<<5)
#define S3C2410_SDICON_BYTEORDER      (1<<4)
#define S3C2410_SDICON_SDIOIRQ        (1<<3)
#define S3C2410_SDICON_RWAITEN        (1<<2)
#define S3C2410_SDICON_FIFORESET      (1<<1)
#define S3C2410_SDICON_CLOCKTYPE      (1<<0)
#define S3C2410_SDICMDCON_ABORT       (1<<12)
#define S3C2410_SDICMDCON_WITHDATA    (1<<11)
#define S3C2410_SDICMDCON_LONGRSP     (1<<10)
#define S3C2410_SDICMDCON_WAITRSP     (1<<9)
#define S3C2410_SDICMDCON_CMDSTART    (1<<8)
#define S3C2410_SDICMDCON_SENDERHOST  (1<<6)
#define S3C2410_SDICMDCON_INDEX       (0x3f)
#define S3C2410_SDICMDSTAT_CRCFAIL    (1<<12)
#define S3C2410_SDICMDSTAT_CMDSENT    (1<<11)
#define S3C2410_SDICMDSTAT_CMDTIMEOUT (1<<10)
#define S3C2410_SDICMDSTAT_RSPFIN     (1<<9)
#define S3C2410_SDICMDSTAT_XFERING    (1<<8)
#define S3C2410_SDICMDSTAT_INDEX      (0xff)
#define S3C2440_SDIDCON_DS_BYTE       (0<<22)
#define S3C2440_SDIDCON_DS_HALFWORD   (1<<22)
#define S3C2440_SDIDCON_DS_WORD       (2<<22)
#define S3C2410_SDIDCON_IRQPERIOD     (1<<21)
#define S3C2410_SDIDCON_TXAFTERRESP   (1<<20)
#define S3C2410_SDIDCON_RXAFTERCMD    (1<<19)
#define S3C2410_SDIDCON_BUSYAFTERCMD  (1<<18)
#define S3C2410_SDIDCON_BLOCKMODE     (1<<17)
#define S3C2410_SDIDCON_WIDEBUS       (1<<16)
#define S3C2410_SDIDCON_DMAEN         (1<<15)
#define S3C2410_SDIDCON_STOP          (1<<14)
#define S3C2440_SDIDCON_DATSTART      (1<<14)
#define S3C2410_SDIDCON_DATMODE                  (3<<12)
#define S3C2410_SDIDCON_BLKNUM        (0x7ff)
/* constants for S3C2410_SDIDCON_DATMODE */
#define S3C2410_SDIDCON_XFER_READY    (0<<12)
#define S3C2410_SDIDCON_XFER_CHKSTART (1<<12)
#define S3C2410_SDIDCON_XFER_RXSTART  (2<<12)
#define S3C2410_SDIDCON_XFER_TXSTART  (3<<12)
#define S3C2410_SDIDCNT_BLKNUM_MASK   (0xFFF)
#define S3C2410_SDIDCNT_BLKNUM_SHIFT  (12)
#define S3C2410_SDIDSTA_RDYWAITREQ    (1<<10)
#define S3C2410_SDIDSTA_SDIOIRQDETECT (1<<9)
#define S3C2410_SDIDSTA_FIFOFAIL      (1<<8)           /* reserved on 2440 */
#define S3C2410_SDIDSTA_CRCFAIL       (1<<7)
#define S3C2410_SDIDSTA_RXCRCFAIL     (1<<6)
#define S3C2410_SDIDSTA_DATATIMEOUT   (1<<5)
#define S3C2410_SDIDSTA_XFERFINISH    (1<<4)
#define S3C2410_SDIDSTA_BUSYFINISH    (1<<3)
#define S3C2410_SDIDSTA_SBITERR       (1<<2)          /* reserved on 2410a/2440 */
#define S3C2410_SDIDSTA_TXDATAON      (1<<1)
#define S3C2410_SDIDSTA_RXDATAON      (1<<0)
#define S3C2440_SDIFSTA_FIFORESET      (1<<16)
#define S3C2440_SDIFSTA_FIFOFAIL
#define S3C2410_SDIFSTA_TFDET
#define S3C2410_SDIFSTA_RFDET
#define S3C2410_SDIFSTA_TFHALF
#define S3C2410_SDIFSTA_TFEMPTY
#define S3C2410_SDIFSTA_RFLAST
#define S3C2410_SDIFSTA_RFFULL
#define S3C2410_SDIFSTA_RFHALF






#define S3C2410_SDIFSTA_COUNTMASK      (0x7f)
#define S3C2410_SDIIMSK_RESPONSECRC    (1<<17)
#define S3C2410_SDIIMSK_CMDSENT        (1<<16)
#define S3C2410_SDIIMSK_CMDTIMEOUT     (1<<15)
#define S3C2410_SDIIMSK_RESPONSEND     (1<<14)
#define S3C2410_SDIIMSK_READWAIT       (1<<13)
#define S3C2410_SDIIMSK_SDIOIRQ        (1<<12)
#define S3C2410_SDIIMSK_FIFOFAIL       (1<<11)
#define S3C2410_SDIIMSK_CRCSTATUS      (1<<10)
#define S3C2410_SDIIMSK_DATACRC        (1<<9)
#define S3C2410_SDIIMSK_DATATIMEOUT    (1<<8)
#define S3C2410_SDIIMSK_DATAFINISH     (1<<7)
#define S3C2410_SDIIMSK_BUSYFINISH     (1<<6)
#define S3C2410_SDIIMSK_SBITERR        (1<<5)          /* reserved 2440/2410a */
#define S3C2410_SDIIMSK_TXFIFOHALF     (1<<4)
#define S3C2410_SDIIMSK_TXFIFOEMPTY    (1<<3)
#define S3C2410_SDIIMSK_RXFIFOLAST     (1<<2)
#define S3C2410_SDIIMSK_RXFIFOFULL     (1<<1)
#define S3C2410_SDIIMSK_RXFIFOHALF     (1<<0)
#endif /* __ASM_ARM_REGS_SDI */