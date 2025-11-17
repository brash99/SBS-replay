#!/usr/bin/env python3

def parse_word(hex_str: str):
    # Allow optional 0x prefix and strip whitespace
    hex_str = hex_str.strip().lower()
    if hex_str.startswith("0x"):
        hex_str = hex_str[2:]

    # Convert to 32-bit integer
    word = int(hex_str, 16) & 0xFFFFFFFF

    # Extract fields (bit 31 is MSB, bit 0 is LSB)
    datatype = (word >> 27) & 0x1F       # bits 31–27 (5 bits)
    group    = (word >> 24) & 0x07       # bits 26–24 (3 bits)
    channel  = (word >> 19) & 0x1F       # bits 23–19 (5 bits)
    edgetype = (word >> 18) & 0x01       # bit 18 (1 bit)
    coarse   = (word >> 8)  & 0x3FF      # bits 17–8 (10 bits)
    twons    = (word >> 7)  & 0x01       # bit 7 (1 bit)
    fine     = word        & 0x7F        # bits 6–0 (7 bits)

    return {
        "datatype": datatype,
        "group": group,
        "channel": channel,
        "edgetype": edgetype,
        "coarse": coarse,
        "twons": twons,
        "fine": fine,
    }

def main():
    hex_word = input("Enter 32-bit hex word (e.g. 0x1234ABCD): ")
    try:
        fields = parse_word(hex_word)
    except ValueError:
        print("Error: invalid hex input.")
        return

    print(f"\nParsed fields from {hex_word}:")
    for name, value in fields.items():
        print(f"  {name:8s} = {value:d}")

if __name__ == "__main__":
    main()
