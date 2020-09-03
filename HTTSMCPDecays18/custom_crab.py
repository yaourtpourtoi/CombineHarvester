from CRABClient.UserUtilities import getUsernameFromCRIC
def custom_crab(config):
    print('>> Customising the crab config')
    config.General.workArea='Impacts-24062020_tol10'
    config.Site.storageSite = 'T2_UK_London_IC'
    config.Site.blacklist=['T2_BE_IIHE']

    config.JobType.maxMemoryMB = 2500

    config.Data.outLFNDirBase='/store/user/{}/{}/'.format(
        getUsernameFromCRIC(), config.General.workArea
    )
    config.JobType.allowUndistributedCMSSW = True
