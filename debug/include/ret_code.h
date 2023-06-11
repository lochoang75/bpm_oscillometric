#ifndef RET_CODE_H
#define RET_CODE_H

typedef enum {
    kRET_OK = 0,
    kRET_INVL,
    kRET_TRY_AGAIN,
    kRET_BUSY,
    kRET_NOMEM,
    kRET_FAILED,
} ret_code_t;
#endif /*RET_CODE_H*/