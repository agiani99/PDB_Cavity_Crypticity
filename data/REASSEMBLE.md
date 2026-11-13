# Reassembling data.zip

This file was split into 11 part(s) for reliable upload. Download all
parts (`data.zip.part001` .. `data.zip.part011`) plus `data.zip.sha256.json`,
then:

## Linux / macOS
```bash
cat data.zip.part* > data.zip
sha256sum data.zip.part*   # compare each against data.zip.sha256.json
```

## Windows
Use `copy /b` (binary concatenation, streams to disk — do not use
`Get-Content`/`Set-Content` for files this large, it loads everything into
memory at once and will exhaust RAM). From PowerShell or cmd.exe:
```
copy /b data.zip.part001+data.zip.part002+data.zip.part003+data.zip.part004+data.zip.part005+data.zip.part006+data.zip.part007+data.zip.part008+data.zip.part009+data.zip.part010+data.zip.part011 data.zip
```
To verify a part's checksum in PowerShell:
```powershell
Get-FileHash data.zip.part001 -Algorithm SHA256
```

## Verifying integrity
Per-part SHA256 checksums are in `data.zip.sha256.json`. Recompute each part's
hash and compare before or after reassembly to confirm nothing was corrupted
in transit.
