import xarray as xr

ds=xr.open_dataset('../runs/short_winter2016_04/DFM_OUTPUT_short_winter2016_04/short_winter2016_04_0000_20151215_000000_his.nc')

## 

# Nice - they include gate discharges automatically
plt.figure(1).clf()
fig,ax=plt.subplots(num=1)
fig.set_size_inches([9,4.8],forward=True)


for gate_idx,gate in enumerate(ds.gategen_name):
    # omit the problem child a14_coyote
    gate=gate.item().decode()
    if gate=='a14_coyote':
        continue
    gate=gate.upper().replace('_','–')
    ax.plot(ds.time, ds.gategen_discharge.isel(gategens=gate_idx),
            label=gate)

ax.legend(fontsize=9)
ax.set_ylabel('Discharge (m3/s)')
ax.axis( (735949.20,735951.81, -19.40, 53.940) )
fig.autofmt_xdate()
ax.axhline(0,color='k',lw=0.8,zorder=-1)


plt.savefig('example-gate-discharges.png')
