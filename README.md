# ESB Sniffer

Enhanced ShockBurst (ESB) is a protocol designed by Nordic. Here we make
a sniffer for ESB packets.

There are 3 types:

1. SB protocol

1. ESB with fixed payload length (`NRF_ESB_PROTOCOL_ESB`);

1. ESB with dynamic payload length (`NRF_ESB_PROTOCOL_ESB_DPL`)

All three protocols are supported. (SB protocol is not tested.)

## Hardware requirements

* [HackRF One](https://greatscottgadgets.com/hackrf/one/)

## Build

Easy, just build `main.cpp` with HackRF library.

### Build on Windows

Since buiding `HackRF` on Windows [looks scaring](https://hackrf.readthedocs.io/en/latest/installing_hackrf_software.html#windows-prerequisites-for-cygwin-mingw-or-visual-studio),
here are detailed steps to build it with MinGW64 which will turn out to be much easier that expected.

* Folder [_libhackrf_](libhackrf) contains source code of the library [itself](https://github.com/greatscottgadgets/hackrf/tree/master/host/libhackrf);

* Folder [_libusb_](libusb) contains header file of [libusb](https://github.com/libusb/libusb/releases/tag/v1.0.27);

    Extract `libusb.dll` for MinGW64 to root directory of this project.

Steps to build:

1. Build `libhackrf`:

    ```sh
    gcc -shared -o libhackrf.dll libhackrf\src\hackrf.c -I libusb libusb-1.0.dll
    ```

1. Build `main.cpp`

    ```sh
    g++ main.cpp libhackrf.dll -I libhackrf\src
    ```

## Command line options

Run `a -h` to check all available options:

```
Basic options:
  -h, --help                      show this help message and exit
BASIC:
  -c, --channel CH                ESB/SB channel (0..100). Default: 2
                                  central frequency = (2040 + ch)MHz
  -f, --freq_hz F                 central frequency in Hz.
                                  k/M/G suffices are supported, such as 2042M, 2.042G.
  --sample_rate R                 sample rate in MHz ({1, 2}). default: 2(MHz)
  -p, --protocol PROTO            protocol (SB/ESB). default: ESB
  --addr_byte_len LEN             address length in bytes ({3, 4, 5}). default: 5
  --crc_byte_len LEN              CRC length in bytes ({1, 2}). default: 2
SB options:
  --payload_len LEN               payload length in bytes ([0..255]). default: 32
ESB options:
  --payload_len LEN               fixed payload length in bytes ([0..255]). default: -1 (dynamic)
  --payload_len_bit_len LEN       payload length field length in bits ({6, 8}). default: 6
RF options:
  --lna_gain GAIN                 LNA (IF) gain in dB ([0, 40], step 8). default: 16
  --vga_gain GAIN                 VGA gain in dB ([0, 62], step 2). default: 16
  +amplifier                      enable RF RX/TX amplifiers. default: disabled
Debug options:
  +verbose                        dump all possible packets.
  --filter AA:BB:...              filter addresses. (default: not filtered)
  --mask   CC:DD:...              mask on address filter. (default: FF:FF:...)
                                  when bit i is given in `filter` and bit i of `mask` is 1,
                                  bit i of address is tested against filter:
                                  if not matched, discard it.
```

## Acknowledgements

* [BTLE](https://github.com/JiaoXianjun/BTLE) project helps a lot;

## References

1. https://infocenter.nordicsemi.com/pdf/nRF24LE1_PS_v1.6.pdf

1. https://devzone.nordicsemi.com/nordic/nordic-blog/b/blog/posts/intro-to-shockburstenhanced-shockburst


