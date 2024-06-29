#include <iostream>
#include <pthread.h>
#include <vector>

#include "hackrf.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <ctype.h>
#include <math.h>

#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <errno.h>
#include <signal.h>

volatile bool is_terminating = false;

#ifdef _WIN32
#include <windows.h>

#ifdef _MSC_VER

#define strtoull _strtoui64
#define snprintf _snprintf

int gettimeofday(struct timeval *tv, void *ignored)
{
    FILETIME ft;
    unsigned __int64 tmp = 0;
    if (NULL != tv)
    {
        GetSystemTimeAsFileTime(&ft);
        tmp |= ft.dwHighDateTime;
        tmp <<= 32;
        tmp |= ft.dwLowDateTime;
        tmp /= 10;
        tmp -= 11644473600000000Ui64;
        tv->tv_sec = (long)(tmp / 1000000UL);
        tv->tv_usec = (long)(tmp % 1000000UL);
    }
    return 0;
}

#endif

#define sleep(a) Sleep((a * 1000))

#endif

#if defined(__GNUC__)
#include <unistd.h>
#include <sys/time.h>
#endif

#ifdef _MSC_VER
BOOL WINAPI sigint_callback_handler(int signum)
{
    if (CTRL_C_EVENT == signum)
    {
        fprintf(stdout, "Caught signal %d\n", signum);
        is_terminating = true;
        return TRUE;
    }
    return FALSE;
}
#else
void sigint_callback_handler(int signum)
{
    if (CTRL_C_EVENT == signum)
    {
        fprintf(stdout, "Caught signal %d\n", signum);
        is_terminating = true;
    }
}
#endif

#define SAMPLE_PER_SYMBOL 4 // 8M sampling rate

#define LEN_BUF_IN_SAMPLE (8 * 4096)
#define LEN_BUF (LEN_BUF_IN_SAMPLE * 2) // must be 2^
#define LEN_BUF_IN_SYMBOL (LEN_BUF_IN_SAMPLE / SAMPLE_PER_SYMBOL)

#define MAX_OVERLAP_BYTE (50)

#define MAX_OVERLAP_SAMPLE (MAX_OVERLAP_BYTE * 8 * SAMPLE_PER_SYMBOL)
#define LEN_BUF_OVERLAP (2 * MAX_OVERLAP_SAMPLE)

typedef int8_t iq_comp_t;

volatile int rx_buf_offset = 0;
volatile iq_comp_t rx_buf[LEN_BUF + LEN_BUF_OVERLAP];

int rx_callback(hackrf_transfer *transfer)
{
    int i;
    int8_t *p = (int8_t *)transfer->buffer;
    for (i = 0; i < transfer->valid_length; i++)
    {
        rx_buf[rx_buf_offset] = p[i];
        rx_buf_offset = (rx_buf_offset + 1) & (LEN_BUF - 1);
        if (rx_buf_offset == LEN_BUF_OVERLAP)
        {
            memcpy((void *)(rx_buf + LEN_BUF), (void *)rx_buf, LEN_BUF_OVERLAP);
        }
    }
    return (0);
}

class RFSource
{
public:
    RFSource()
        : device(nullptr),
          freq_hz(1000000),
          sample_rate(1000000),
          vga_gain(16),
          lna_gain(16),
          amp_enable(false)
    {
    }

    virtual ~RFSource()
    {
        close();
    }

public:
    int run(hackrf_sample_block_cb_fn callback)
    {
        int result;

        result = init();
        if (result != 0)
            return result;

        result = open();
        if (result != 0)
            return result;

        result = hackrf_start_rx(device, callback, NULL);
        if (result != HACKRF_SUCCESS)
        {
            printf("RFSource::run(): hackrf_start_rx() failed: %s (%d)\n", hackrf_error_name((hackrf_error)result), result);
            return -1;
        }

        return 0;
    }

protected:
    int init()
    {
        int result = hackrf_init();
        if (result != HACKRF_SUCCESS)
        {
            printf("RFSource::init(): hackrf_init() failed: %s (%d)\n", hackrf_error_name((hackrf_error)result), result);
            return -1;
        }

        return 0;
    }

    int open()
    {
        int result;

        result = hackrf_open(&device);
        if (result != HACKRF_SUCCESS)
        {
            printf("RFSource::open(): hackrf_open() failed: %s (%d)\n", hackrf_error_name((hackrf_error)result), result);
            return -1;
        }

        printf("Setting central frequency to %zdHz\n", freq_hz);
        result = hackrf_set_freq(device, freq_hz);
        if (result != HACKRF_SUCCESS)
        {
            printf("RFSource::open(): hackrf_set_freq() failed: %s (%d)\n", hackrf_error_name((hackrf_error)result), result);
            return -1;
        }

        result = hackrf_set_sample_rate(device, (double)sample_rate);
        if (result != HACKRF_SUCCESS)
        {
            printf("RFSource::open(): hackrf_set_sample_rate() failed: %s (%d)\n", hackrf_error_name((hackrf_error)result), result);
            return -1;
        }

        result = hackrf_set_baseband_filter_bandwidth(device, sample_rate / 2);
        if (result != HACKRF_SUCCESS)
        {
            printf("RFSource::open(): hackrf_set_baseband_filter_bandwidth() failed: %s (%d)\n", hackrf_error_name((hackrf_error)result), result);
            return -1;
        }

        printf("Setting VGA gain to %d\n", vga_gain);
        result = hackrf_set_vga_gain(device, vga_gain);
        if (result != HACKRF_SUCCESS)
        {
            printf("RFSource::open(): hackrf_set_vga_gain() failed: %s (%d)\n", hackrf_error_name((hackrf_error)result), result);
            return -1;
        }

        printf("Setting LNA gain to %d\n", lna_gain);
        result = hackrf_set_lna_gain(device, lna_gain);
        if (result != HACKRF_SUCCESS)
        {
            printf("RFSource::open(): hackrf_set_lna_gain() failed: %s (%d)\n", hackrf_error_name((hackrf_error)result), result);
            return -1;
        }

        printf(amp_enable ? "Enabling amp\n" : "Disabling amp\n");
        result = hackrf_set_amp_enable(device, amp_enable ? 1 : 0);
        if (result != HACKRF_SUCCESS)
        {
            printf("RFSource::open(): hackrf_set_amp_enable() failed: %s (%d)\n", hackrf_error_name((hackrf_error)result), result);
            return -1;
        }

        return 0;
    }

    int close()
    {
        int result;

        if (device == nullptr)
            return -1;

        result = hackrf_stop_rx(device);
        if (result != HACKRF_SUCCESS)
        {
            printf("close_board: hackrf_stop_rx() failed: %s (%d)\n", hackrf_error_name((hackrf_error)result), result);
            return -1;
        }

        result = hackrf_close(device);
        if (result != HACKRF_SUCCESS)
        {
            printf("close_board: hackrf_close() failed: %s (%d)\n", hackrf_error_name((hackrf_error)result), result);
            return -1;
        }

        hackrf_exit();
        printf("hackrf_exit() done\n");

        return 0;
    }

public:
    hackrf_device *device;
    uint64_t freq_hz;
    uint64_t sample_rate;
    int vga_gain;
    int lna_gain;
    bool amp_enable;
};

void save_phy_sample(iq_comp_t *iq_sample, int num_iq_sample, char *filename)
{
    int i;

    FILE *fp = fopen(filename, "w");
    if (fp == NULL)
    {
        printf("save_phy_sample: fopen failed!\n");
        return;
    }

    for (i = 0; i < num_iq_sample; i++)
    {
        if (i % 64 == 0)
        {
            fprintf(fp, "\n");
        }
        fprintf(fp, "%d, ", iq_sample[i]);
    }
    fprintf(fp, "\n");

    fclose(fp);
}

void dem_bits(const iq_comp_t *rxp, int num_bits, uint32_t &out_packed_bits)
{
    int i, j;
    int I0, Q0, I1, Q1;
    int sample_idx = 0;

    uint32_t acc = 0;

    for (i = 0; i < num_bits; i++)
    {
        I0 = rxp[sample_idx];
        Q0 = rxp[sample_idx + 1];
        I1 = rxp[sample_idx + 2];
        Q1 = rxp[sample_idx + 3];

        uint8_t bit_decision = (I0 * Q1 - I1 * Q0) > 0 ? 1 : 0;

        acc <<= 1;
        acc |= bit_decision;

        sample_idx = sample_idx + SAMPLE_PER_SYMBOL * 2;
    }

    out_packed_bits = acc;
}

void dem_byte(const iq_comp_t *rxp, int num_byte, uint8_t *out_byte)
{
    int i, j;
    int I0, Q0, I1, Q1;
    int sample_idx = 0;

    for (i = 0; i < num_byte; i++)
    {
        out_byte[i] = 0;
        for (j = 0; j < 8; j++)
        {
            I0 = rxp[sample_idx];
            Q0 = rxp[sample_idx + 1];
            I1 = rxp[sample_idx + 2];
            Q1 = rxp[sample_idx + 3];

            uint8_t bit_decision = (I0 * Q1 - I1 * Q0) > 0 ? 1 : 0;
            out_byte[i] = out_byte[i] | (bit_decision << (7 - j));

            sample_idx = sample_idx + SAMPLE_PER_SYMBOL * 2;
        }
    }
}

void dem_byte_le(const iq_comp_t *rxp, int num_byte, uint8_t *out_byte)
{
    int i, j;
    int I0, Q0, I1, Q1;
    int sample_idx = 0;

    for (i = 0; i < num_byte; i++)
    {
        out_byte[i] = 0;
        for (j = 0; j < 8; j++)
        {
            I0 = rxp[sample_idx];
            Q0 = rxp[sample_idx + 1];
            I1 = rxp[sample_idx + 2];
            Q1 = rxp[sample_idx + 3];

            uint8_t bit_decision = (I0 * Q1 - I1 * Q0) > 0 ? 1 : 0;
            out_byte[i] = out_byte[i] | (bit_decision << j);

            sample_idx = sample_idx + SAMPLE_PER_SYMBOL * 2;
        }
    }
}

template <int num_sync_bits, int sample_per_symbol>
class PhaseSync
{
public:
    PhaseSync(uint32_t sync, const uint8_t *sync_mask_bytes = nullptr)
    {
        init(sync, sync_mask_bytes);
    }

    int sync(const iq_comp_t *rxp, const int search_len_in_symbol)
    {
        int phase_idx;
        const int dem_buf_len = num_sync_bits;
        int dem_buf_offset = 0;
        memset(dem_buf, 0, sizeof(dem_buf));

        for (int i = 0; i < search_len_in_symbol * sample_per_symbol * 2; i += sample_per_symbol * 2)
        {
            int sp = ((dem_buf_offset - dem_buf_len + 1) & (dem_buf_len - 1));

            for (int j = 0; j < sample_per_symbol * 2; j += 2)
            {
                int i0 = rxp[i + j];
                int q0 = rxp[i + j + 1];
                int i1 = rxp[i + j + 2];
                int q1 = rxp[i + j + 3];

                phase_idx = j / 2;
                dem_buf[phase_idx][dem_buf_offset] = (i0 * q1 - i1 * q0) > 0 ? 1 : 0;

                int k = sp;
                bool mismatch = false;
                for (int p = 0; p < dem_buf_len; p++)
                {
                    if (sync_bits_mask[p] && (dem_buf[phase_idx][k] != sync_bits[p]))
                    {
                        mismatch = true;
                        break;
                    }
                    k = (k + 1) & (dem_buf_len - 1);
                }

                if (!mismatch)
                {
                    return i + j - (dem_buf_len - 1) * sample_per_symbol * 2;
                }
            }

            dem_buf_offset = (dem_buf_offset + 1) & (dem_buf_len - 1);
        }

        return -1;
    }

    int get_num_sync_bits(void)
    {
        return num_sync_bits;
    }

public:
    uint8_t sync_bits[num_sync_bits];
    uint8_t sync_bits_mask[num_sync_bits];
    uint8_t dem_buf[sample_per_symbol][num_sync_bits * 2];

private:
    void init(uint32_t sync, const uint8_t *sync_mask_bytes)
    {
        int i = 0;
        for (int i  = 0; i < num_sync_bits; i++)
        {
            uint32_t b = sync >> (num_sync_bits - 1 - i);
            uint8_t m = sync_mask_bytes ? *sync_mask_bytes++ : 0xff;
            this->sync_bits[i] = b & 1;
            this->sync_bits_mask[i] = m & 1;
        }
    }
};

void disp_hex(uint8_t *hex, int num_hex, bool nl = true)
{
    int i;
    for (i = 0; i < num_hex; i++)
    {
        printf("%02x", hex[i]);
    }
    if (nl)
        printf("\n");
}

template <class crc_t, crc_t POLYNOMIAL> class CRC
{
public:
    CRC():
        WIDTH(8 * sizeof(crc_t)),
        TOPBIT(1 << ((8 * sizeof(crc_t)) - 1)),
        remainder(0)
    {
        crc_t  remainder;

        for (int dividend = 0; dividend < 256; ++dividend)
        {
            remainder = dividend << (WIDTH - 8);

            for (uint8_t bit = 8; bit > 0; --bit)
            {
                if (remainder & TOPBIT)
                {
                    remainder = (remainder << 1) ^ POLYNOMIAL;
                }
                else
                {
                    remainder = (remainder << 1);
                }
            }

            crc_table[dividend] = remainder;
        }
    }

    void restart(crc_t init = 0)
    {
        remainder = init;
    }

    crc_t update(const uint8_t *message, int byte_len)
    {
        uint8_t data;

        for (int byte = 0; byte < byte_len; ++byte)
        {
            data = message[byte] ^ (remainder >> (WIDTH - 8));
            remainder = crc_table[data] ^ (remainder << 8);
        }

        return remainder;
    }

    // starting from lowest bit
    crc_t update(uint32_t value, int bit_len)
    {
        for (int i = 0; i < bit_len; i++)
        {
            uint8_t bit = (value >> (bit_len - 1 - i)) & 1; // (value >> (bit_len - 1 - i)) & 1;// (value >> (i)) & 1;
            remainder ^= ((crc_t)bit) << (WIDTH - 1);

            if (remainder & TOPBIT)
            {
                remainder <<= 1;
                remainder ^= POLYNOMIAL;
            }
            else
            {
                remainder <<= 1;
            }
        }
        return remainder;
    }

    crc_t get_result(void)
    {
        return remainder;
    }
protected:
    const int WIDTH;
    const int TOPBIT;
    crc_t crc_table[256];
    crc_t remainder;
};

static uint8_t crc_8(uint8_t *addr, uint32_t num, uint8_t init, uint16_t poly)
{
    uint8_t data;
    uint8_t crc = init;
    int i;
    for (; num > 0; num--)
    {
        data = *addr++;

        crc = crc ^ data;
        for (i = 0; i < 8; i++)
        {
            if (crc & 0x80)
            {
                crc = (crc << 1) ^ poly;
            }
            else
            {
                crc <<= 1;
            }
        }
    }
    return crc;
}

#define MHz 1000000ul

class ESBDecoder
{
public:
    ESBDecoder(int addr_byte_len = 5, int payload_len_bit_len = 6, int crc_byte_len = 2, bool verbose = false)
        : sync55(0x55), syncAA(0xAA),
          addr_byte_len(addr_byte_len),
          payload_len_bit_len(payload_len_bit_len),
          crc_byte_len(crc_byte_len),
          verbose(verbose)
    {
    }

    void decode(const iq_comp_t *rxp_in, int buf_len)
    {
        const iq_comp_t *rxp = rxp_in;

        const int total_symbol = buf_len / (SAMPLE_PER_SYMBOL * 2); // 2 for IQ
        int num_symbol_left = total_symbol;

        while (num_symbol_left > MAX_OVERLAP_SAMPLE)
        {
            int hit_idx = sync(rxp, num_symbol_left);
            if (hit_idx == -1)
                break;

            rxp += hit_idx + sync55.get_num_sync_bits() * 2 * SAMPLE_PER_SYMBOL; // move to beginning of PDU header
            const iq_comp_t *pkt_start = rxp;

            dem_byte(rxp, addr_byte_len + 2, addr);
            rxp += addr_byte_len * 8 * SAMPLE_PER_SYMBOL * 2;

            uint32_t acc = *(uint32_t *)addr;

            uint32_t payload_len = 0;
            uint32_t pid = 0;
            uint32_t no_ack = 0;

            dem_bits(rxp, payload_len_bit_len, payload_len);
            rxp += payload_len_bit_len * SAMPLE_PER_SYMBOL * 2;

            dem_bits(rxp, 2, pid);
            rxp += 2 * SAMPLE_PER_SYMBOL * 2;

            dem_bits(rxp, 1, no_ack);
            rxp += 1 * SAMPLE_PER_SYMBOL * 2;

            uint8_t len = payload_len;

            uint8_t crc_value[2] = {0};
            dem_byte(rxp, len, payload);
            rxp += len * 8 * SAMPLE_PER_SYMBOL * 2;
            dem_byte(rxp, crc_byte_len, crc_value);
            rxp += crc_byte_len * 8 * SAMPLE_PER_SYMBOL * 2;

            uint16_t air_crc = crc_byte_len == 2 ? (crc_value[0] << 8) | crc_value[1] : crc_value[0];
            uint16_t cal_crc = calc_crc(addr, len, (uint8_t)pid, (uint8_t)no_ack, payload);
            if (air_crc == cal_crc)
            {
                disp_hex(addr, addr_byte_len, false);
                if (len > 0)
                {
                    printf(": PID=%d NO_ACK=%d LEN=%03d DATA=", pid, no_ack, len);
                    disp_hex(payload, len);
                }
                else
                    printf(": PID=%d NO_ACK=%d ACK\n", pid, no_ack);
            }
            else
            {
                if (verbose)
                {
                    printf("BAD: ");
                    disp_hex(addr, addr_byte_len, false);
                    printf(": PID=%d NO_ACK=%d LEN=%03d DATA=", pid, no_ack, len);
                    disp_hex(payload, len, false);
                    printf(" CRC %04x <> %04x\n", air_crc, cal_crc);
                }
                rxp = pkt_start;
            }

            num_symbol_left = total_symbol - (rxp - rxp_in) / (SAMPLE_PER_SYMBOL * 2);
        }
    }

protected:
    PhaseSync<8, SAMPLE_PER_SYMBOL> sync55;
    PhaseSync<8, SAMPLE_PER_SYMBOL> syncAA;
    const int addr_byte_len;       // [3,5]
    const int payload_len_bit_len; // {6, 8}
    const int crc_byte_len;        // {0, 1, 2}
    const bool verbose;

    uint8_t addr[10];
    uint8_t payload[300]; // payload and crc

    CRC<uint8_t, 0x7> crc8;
    CRC<uint16_t, 0x1021> crc16;

protected:

    int sync(PhaseSync<8, SAMPLE_PER_SYMBOL> &syncer, const iq_comp_t *rxp, const int search_len_in_symbol, uint32_t expected_addr_bit)
    {
        int hit = syncer.sync(rxp, search_len_in_symbol);
        if (hit == -1)
            return hit;

        rxp += hit + syncer.get_num_sync_bits() * 2 * SAMPLE_PER_SYMBOL; // move to beginning of PDU header

        uint32_t addr;
        dem_bits(rxp, 1, addr);

        return addr == expected_addr_bit ? hit : -1;
    }

    int sync(const iq_comp_t *rxp, const int search_len_in_symbol)
    {
        int hit55 = sync(sync55, rxp, search_len_in_symbol, 0);
        int hitAA = sync(syncAA, rxp, search_len_in_symbol, 1);

        if (hit55 == -1) return hitAA;
        if (hitAA == -1) return hit55;

        return hit55 < hitAA ? hit55 : hitAA;
    }

    uint16_t calc_crc(const uint8_t *addr, uint8_t len, uint8_t pid, uint8_t no_ack, const uint8_t *payload)
    {
        if (crc_byte_len == 2)
        {
            auto &crc = crc16;
            crc.restart(0xffff);
            crc.update(addr, addr_byte_len);
            crc.update(len, payload_len_bit_len);
            crc.update(pid, 2);
            crc.update(no_ack, 1);
            crc.update(payload, len);
            return crc.get_result();
        }
        else
        {
            auto &crc = crc8;
            crc.restart(0xff);
            crc.update(addr, addr_byte_len);
            crc.update(len, payload_len_bit_len);
            crc.update(pid, 2);
            crc.update(no_ack, 1);
            crc.update(payload, len);
            return crc.get_result();
        }
    }
};

struct Args
{
    uint64_t freq_hz = 2402 * MHz;
    int addr_byte_len = 5;
    int payload_len_bit_len = 6;
    int crc_byte_len = 2;
    int lna_gain = 16;
    int vga_gain = 16;
    int sample_rate = 2;
    bool amp_enable = false;
    bool verbose = false;
    bool show_help = false;
};

void usage(const std::string &prog)
{
    std::cout << "Usage: " << prog << " [options]\n"
              << "\n"
              << "Basic options:\n"
              << "  -h, --help                      show this help message and exit\n"
              << "BASIC:\n"
              << "  -f, --freq_hz F                 frequency in Hz. default: 2042M\n"
              << "                                  k/M/G suffices are supported, such as 2042M, 2.042G.\n"
              << "  --sample_rate R                 sample rate in MHz ({1, 2}). default: 2(MHz)\n"
              << "  --addr_byte_len LEN             address length in bytes ({3, 4, 5}). default: 5\n"
              << "  --payload_len_bit_len LEN       payload length field length in bits ({6, 8}). default: 6\n"
              << "  --crc_byte_len LEN              CRC length in bytes ({1, 2}). default: 2\n"
              << "RF options:\n"
              << "  --lna_gain GAIN                 LNA (IF) gain in dB ([0, 40], step 8). default: 16\n"
              << "  --vga_gain GAIN                 VGA gain in dB ([0, 62], step 2). default: 16\n"
              << "  +amplifier                      enable RF RX/TX amplifiers. default: disabled\n"
              << "  +verbose                        dump all possible packets.\n";
}

static uint64_t parse_freq(const char *s)
{
    double r = atof(s);
    char l = s[strlen(s) - 1];
    switch (l)
    {
    case 'k':
        r *= 1000;
        break;
    case 'M':
        r *= 1000000;
        break;
    case 'G':
        r *= 1000000000;
        break;
    default:
        break;
    }
    return (uint64_t)r;
}

static size_t parse_args(Args &args, const std::vector<std::string> &argv)
{
    const size_t argc = argv.size();

    #define handle_para0(fmt1, field, f)        \
        else if ((strcmp(arg, fmt1) == 0))      \
        {                                                                   \
            c++;                                                            \
            if (c < argc)                                                   \
                args.field = f(argv[c].c_str());                            \
        }

    #define handle_param(fmt1, fmt2, field, f)    \
        else if ((strcmp(arg, fmt1) == 0) || (strcmp(arg, fmt2) == 0))      \
        {                                                                   \
            c++;                                                            \
            if (c < argc)                                                   \
                args.field = f(argv[c].c_str());                            \
        }

    size_t c = 1;

    try
    {
        while (c < argc)
        {
            const char *arg = argv[c].c_str();
            if ((strcmp(arg, "--help") == 0) || (strcmp(arg, "-h") == 0))
            {
                args.show_help = true;
            }
            else if (strcmp(arg, "+verbose") == 0)
            {
                args.verbose = true;
            }
            else if (strcmp(arg, "+amplifier") == 0)
            {
                args.amp_enable = true;
            }
            handle_param("--freq_hz",               "-f", freq_hz,              parse_freq)
            handle_para0("--addr_byte_len",               addr_byte_len,        std::stoi)
            handle_para0("--payload_len_bit_len",         payload_len_bit_len,  std::stoi)
            handle_para0("--crc_byte_len",                crc_byte_len,         std::stoi)
            handle_para0("--sample_rate",                 sample_rate,          std::stoi)
            handle_para0("--lna_gain",                    lna_gain,             std::stoi)
            handle_para0("--vga_gain",                    vga_gain,             std::stoi)
            else
                break;

            c++;
        }
    }
    catch(const std::exception& e)
    {
        std::cerr << e.what() << '\n';
        return c;
    }

    return c;
}

int main(int argc, char **argv)
{
    int phase = 0;

    Args args;

    std::vector<std::string> all_args;
    for (int i = 0; i < argc; i++)
        all_args.push_back(argv[i]);

    const int count = parse_args(args, all_args);

    args.lna_gain = (args.lna_gain + 4) / 8 * 8;
    args.vga_gain = (args.vga_gain + 1) / 2 * 2;

    if (args.show_help
        || (args.sample_rate < 1) || (args.sample_rate > 2)
        || (args.addr_byte_len < 3) || (args.addr_byte_len > 5)
        || (args.crc_byte_len < 1) || (args.crc_byte_len > 2)
        || (args.payload_len_bit_len < 2)|| (args.payload_len_bit_len > 8))
    {
        usage(all_args[0]);
        return 0;
    }

    if (count < all_args.size())
    {
        std::cerr << "Unknown arguments:";
        for (auto i = count; i < all_args.size(); i++)
        {
            std::cerr << " " << all_args[i];
        }
        std::cerr << std::endl;

        exit(EXIT_FAILURE);
    }

#ifdef _MSC_VER
    SetConsoleCtrlHandler((PHANDLER_ROUTINE)sighandler, TRUE);
#else
    signal(SIGINT, &sigint_callback_handler);
    signal(SIGILL, &sigint_callback_handler);
    signal(SIGFPE, &sigint_callback_handler);
    signal(SIGSEGV, &sigint_callback_handler);
    signal(SIGTERM, &sigint_callback_handler);
    signal(SIGABRT, &sigint_callback_handler);
#endif

    RFSource source;
    source.freq_hz = args.freq_hz;
    source.sample_rate = SAMPLE_PER_SYMBOL * args.sample_rate * MHz;
    source.lna_gain = args.lna_gain;
    source.vga_gain = args.vga_gain;
    source.amp_enable = args.amp_enable;

    ESBDecoder *decoder = new ESBDecoder(args.addr_byte_len, args.payload_len_bit_len, args.crc_byte_len, args.verbose);

    if (source.run(rx_callback) != 0)
        return -1;

    while (!is_terminating)
    {
        iq_comp_t *rxp = nullptr;

        int offset = rx_buf_offset - LEN_BUF_OVERLAP;

        if ((offset >= 0) && (rx_buf_offset < LEN_BUF / 2) && (phase == 1))
        {
            phase = 0;
            rxp = (iq_comp_t *)rx_buf + (LEN_BUF / 2);
        }

        if ((offset >= LEN_BUF / 2) && (phase == 0))
        {
            phase = 1;
            rxp = (iq_comp_t *)rx_buf;
        }

        if (rxp)
        {
            decoder->decode(rxp, LEN_BUF / 2 + MAX_OVERLAP_SAMPLE);
            fflush(stdout);
        }
    }

    printf("Exit main loop ...\n");

    return 0;
}
