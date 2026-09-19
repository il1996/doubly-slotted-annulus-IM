# Solve and export the four remaining voltages, one COM session per step.
$T = "<home>\Desktop\claude\T2_ansys"
Set-Location $T
Add-Content "$T\batch_status.txt" "COM CHAIN START $(Get-Date -Format 'HH:mm:ss')"
foreach ($d in @("NL_V550", "NL_V410", "NL_V275", "NL_V140")) {
    if (Test-Path "$T\exports\${d}_NL_wave.tab") { Add-Content "$T\batch_status.txt" "$d already exported, skipped"; continue }
    Add-Content "$T\batch_status.txt" "$d solve begin $(Get-Date -Format 'HH:mm:ss')"
    Remove-Item "$T\IM_18kW_690V_noload_sweep.aedt.lock" -ErrorAction SilentlyContinue
    python "$T\solve_one_com.py" $d 2>&1 | Out-File "$T\solve_$d.txt" -Encoding utf8
    Add-Content "$T\batch_status.txt" "$d solve end $(Get-Date -Format 'HH:mm:ss')"
    Start-Sleep -Seconds 10
    Remove-Item "$T\IM_18kW_690V_noload_sweep.aedt.lock" -ErrorAction SilentlyContinue
    python "$T\noload_export_com.py" $d --close 2>&1 | Out-File "$T\export_$d.txt" -Encoding utf8
    Add-Content "$T\batch_status.txt" "$d export end $(Get-Date -Format 'HH:mm:ss') exported=$(Test-Path "$T\exports\${d}_NL_wave.tab")"
}
Add-Content "$T\batch_status.txt" "COM CHAIN DONE $(Get-Date -Format 'HH:mm:ss')"
