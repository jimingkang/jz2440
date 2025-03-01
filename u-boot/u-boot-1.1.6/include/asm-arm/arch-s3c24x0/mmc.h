#ifndef __MMC_PXA_P_H     
#define __MMC_PXA_P_H     

#include <asm/arch/regs-sdi.h>	
#define MMC_DEFAULT_RCA	(1<<16)
#define MMC_BLOCK_SIZE	512
#define MMC_CMD_RESET	0
#define MMC_CMD_SEND_OP_COND	1
#define MMC_CMD_ALL_SEND_CID	2
#define MMC_CMD_SET_RCA	3
#define MMC_CMD_SELECT_CARD	7
#define MMC_CMD_IF_COND                 8
#define MMC_CMD_S ND_CSD	9
#define MMC_CMD_SEND_CID	10
#define MMC_CMD_SEND_STATUS	13
#define MMC_CMD_SET_BLOCKLEN	16
#define MMC_CMD_READ_BLOCK	17
#define MMC_CMD_RD_BLK_MULTI	18
#define MMC_CMD_WRITE_BLOCK	24
#define MMC_MAX_BLOCK_SIZE	512
#define MMC_R1_IDLE_STATE              0x01	
#define MMC_R1_ERASE_STATE	0x02
#define MMC_R1_ILLEGAL_CMD	0x04
#define MMC_R1_COM_CRC_ERR	0x08
#define MMC_R1_ERASE_SEQ_ERR	0x01
#define MMC_R1_ADDR_ERR	0x02
#define MMC_R1_PARAM_ERR	0x04
#define MMC_R1B_WP_ERASE_SKIP	0x0002
#define MMC_R1B_ERR	0x0004
#define MMC_R1B_CC_ERR	0x0008
#define MMC_R1B_CARD_ECC_ERR	0x0010
#define MMC_R1B_WP_VIOLATION	0x0020
#define MMC_R1B_ERASE_PARAM	0x0040
#define MMC_R1B_OOR	0x0080
#define MMC_R1B_IDLE_STATE	0x0100
#define MMC_R1B_ERASE_RESET	0x0200
#define MMC_R1B_ILLEGAL_CMD	0x0400
#define MMC_R1B_COM_CRC_ERR	0x0800
#define MMC_R1B_ERASE_SEQ_ERR	0x1000
#define MMC_R1B_ADDR_ERR	0x2000
#define MMC_R1B_PARAM_ERR	0x4000
typedef struct mmc_cid
{
/* FIXME: BYTE_ORDER */
uchar  year:4,
month:4; 
uchar sn[3]; uchar fwrev:4,  hwrev:4;
uchar  name[6];
	uchar id[3];
     } mmc_cid_t;

typedef struct mmc_csd
{
uchar    ecc:2,file_format:2,tmp_write_protect:1,perm_write_protect:1,copy:1,file_format_grp:1;
uint64_t content_prot_app:1,rsvd3:4,write_bl_partial:1,write_bl_len:4,r2w_factor:3,default_ecc:2,wp_grp_enable:1,wp_grp_size:5,erase_grp_mult:5,erase_grp_size:5,c_size_mult1:3,vdd_w_curr_max:3,vdd_w_curr_min:3,vdd_r_curr_max:3,vdd_r_curr_min:3,
c_size:12,
rsvd2:2,
dsr_imp:1,
read_blk_misalign:1,
write_blk_misalign:1,
read_bl_partial:1,
read_bl_len:4,
ccc:12;
uchar  tran_speed;
uchar  nsac;
uchar  taac;
uchar rsvd1:2,
  spec_vers:4,
 csd_structure:2;

} mmc_csd_t;	
#endif /* __MMC_PXA  P  H     */